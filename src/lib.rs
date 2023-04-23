use std::{fmt};
use std::io::{Read, Cursor};
use byteorder_pack::UnpackFrom;

type Result<T> = std::result::Result<T, DmapError>;

#[derive(Debug, Clone)]
struct DmapError {
    details: String
}

impl DmapError {
    fn new(msg: String) -> DmapError {
        DmapError{details: msg.to_string()}
    }
}

impl fmt::Display for DmapError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}", self.details)
    }
}

#[derive(Debug)]
#[repr(C)]
pub enum DmapType {
    DMAP,
    CHAR(i8),
    SHORT(i16),
    INT(i32),
    FLOAT(f32),
    DOUBLE(f64),
    STRING(String),
    LONG(i64),
    UCHAR(u8),
    USHORT(u16),
    UINT(u32),
    ULONG(u64),
}

impl DmapType {
    fn all_keys() -> Vec<i8> { vec![0, 1, 2, 3, 4, 8, 9, 10, 16, 17, 18, 19] }

    fn get_num_bytes(self) -> u64 {
        match self {
            DmapType::CHAR {..}   => 1,
            DmapType::SHORT {..}  => 2,
            DmapType::INT {..}    => 4,
            DmapType::FLOAT {..}  => 4,
            DmapType::DOUBLE {..} => 8,
            DmapType::LONG {..}   => 8,
            DmapType::UCHAR {..}  => 1,
            DmapType::USHORT {..} => 2,
            DmapType::UINT {..}   => 4,
            DmapType::ULONG {..}  => 8,
            _                     => 0
        }
    }

    fn get_type_from_fmt(fmt: char) -> Result<DmapType> {
        match fmt {
            'c' => Ok(DmapType::CHAR(0)),
            'h' => Ok(DmapType::SHORT(0)),
            'i' => Ok(DmapType::INT(0)),
            'f' => Ok(DmapType::FLOAT(0.0)),
            'd' => Ok(DmapType::DOUBLE(0.0)),
            's' => Ok(DmapType::STRING("".to_string())),
            'q' => Ok(DmapType::LONG(0)),
            'B' => Ok(DmapType::UCHAR(0)),
            'H' => Ok(DmapType::USHORT(0)),
            'I' => Ok(DmapType::UINT(0)),
            'Q' => Ok(DmapType::ULONG(0)),
            _   => Err(DmapError::new(format!("Invalid format for DMAP type: {}", fmt)))
        }
    }

    fn get_type_from_key(key: i8) -> Result<DmapType> {
        match key {
            0  => Ok(DmapType::DMAP),
            1  => Ok(DmapType::CHAR(0)),
            2  => Ok(DmapType::SHORT(0)),
            3  => Ok(DmapType::INT(0)),
            4  => Ok(DmapType::FLOAT(0.0)),
            8  => Ok(DmapType::DOUBLE(0.0)),
            9  => Ok(DmapType::STRING("".to_string())),
            10 => Ok(DmapType::LONG(0)),
            16 => Ok(DmapType::UCHAR(0)),
            17 => Ok(DmapType::USHORT(0)),
            18 => Ok(DmapType::UINT(0)),
            19 => Ok(DmapType::ULONG(0)),
            _  => Err(DmapError::new(format!("Invalid key for DMAP type: {}", key)))
        }
    }
}

struct RawDmapScalar {
    data: DmapType,
    name: String,
    mode: i8,
}

struct RawDmapArray {
    dmap_type: DmapType,
    name: Box<str>,
    mode: Box<str>,
    dimension: u32,
    arr_dimensions: Vec<u32>,
    data: Vec<u8>,
    data_type_fmt: char
}

struct RawDmapRecord {
    num_scalars: i32,
    num_arrays: i32,
    scalars: Vec<RawDmapScalar>,
    arrays: Vec<RawDmapArray>
}

struct RawDmapRead {
    cursor: Cursor<Vec<u8>>,
    dmap_records: Vec<RawDmapRecord>,
}

impl RawDmapRead {

    fn new(dmap_data: &mut &impl Read) -> Result<RawDmapRead> {
        let mut buffer: Vec<u8> = vec![];

        dmap_data.read_to_end(&mut buffer)
            .map_err(|_| DmapError::new("Could not read data".to_string()))?;

        let mut cursor = Cursor::new(buffer);
        let mut dmap_read = RawDmapRead{cursor, dmap_records: vec![]};

        // TODO: Test initial data integrity

        while dmap_read.cursor.position() < dmap_read.cursor.get_ref().len() as u64 {
            dmap_read.parse_record()?;
        }
        Ok(dmap_read)
    }

    fn parse_record(&mut self) -> Result<()> {
        let bytes_already_read = self.cursor.position();
        let code = match self.read_data(DmapType::INT(0))? {
            DmapType::INT(i) => Ok(i),
            _ => Err(DmapError::new("PARSE RECORD: Invalid code".to_string()))
        }?;
        let size = match self.read_data(DmapType::INT(0))? {
            DmapType::INT(i) => Ok(i),
            _ => Err(DmapError::new("PARSE RECORD: Invalid size".to_string()))
        }?;

        // adding 8 bytes because code and size are part of the record.
        if size as u64 > self.cursor.get_ref().len() as u64
            - self.cursor.position()
            + 2 * DmapType::INT(0).get_num_bytes() {
            return Err(DmapError::new("PARSE RECORD: Integrity check shows record size bigger than \
                remaining buffer. Data is likely corrupted".to_string()))
        }
        else if size <= 0 {
            return Err(DmapError::new("PARSE RECORD: Integrity check shows record size <= 0. \
                Data is likely corrupted".to_string()))
        }

        let num_scalars = match self.read_data(DmapType::INT(0))? {
            DmapType::INT(i) => Ok(i),
            _ => Err(DmapError::new("PARSE RECORD: Invalid number of scalars".to_string()))
        }?;
        let num_arrays = match self.read_data(DmapType::INT(0))? {
            DmapType::INT(i) => Ok(i),
            _ => Err(DmapError::new("PARSE RECORD: Invalid number of arrays".to_string()))
        }?;

        if num_scalars <= 0 {
            return Err(DmapError::new("PARSE RECORD: Number of scalers is 0 or negative.".to_string()))
        }
        else if num_arrays <= 0 {
            return Err(DmapError::new("PARSE RECORD: Number of arrays is 0 or negative.".to_string()))
        }
        else if num_scalars + num_arrays > size {
            return Err(DmapError::new("PARSE RECORD: Invalid number of record elements. \
                Array or scaler field is likely corrupted.".to_string()))
        };

        let mut scalars: Vec<RawDmapScalar> = vec![];
        for n in 0..num_scalars {
            scalars.push(self.parse_scalar()?);
        }
        let mut arrays: Vec<RawDmapArray> = vec![];
        for i in 0..num_arrays {
            arrays.push(self.parse_array(size));
        }

        if self.cursor.position() - bytes_already_read != size as u64 {
            return Err(DmapError::new(format!(
                "PARSE RECORD: Bytes read {} does not match the records size field {}",
                self.cursor.position() - bytes_already_read, size)))
        }

        self.dmap_records
            .push(RawDmapRecord{num_scalars, scalars, num_arrays, arrays});
        Ok(())
    }

    fn parse_scalar(&mut self) -> Result<RawDmapScalar> {
        let mode = 6;
        let name = match self.read_data(DmapType::STRING("".to_string()))? {
            DmapType::STRING(s) => Ok(s),
            _ => Err(DmapError::new("PARSE SCALAR: Invalid scalar name".to_string()))
        }?;
        let data_type_key = match self.read_data(DmapType::CHAR(0))? {
            DmapType::CHAR(c) => Ok(c),
            _ => Err(DmapError::new("PARSE SCALAR: Invalid data type".to_string()))
        }?;

        if !DmapType::all_keys().contains(&data_type_key) {
            return Err(DmapError::new("PARSE SCALAR: Data type is corrupted. Record is likely corrupted".to_string()))
        }

        let data_type = DmapType::get_type_from_key(data_type_key)?;

        let data = match data_type {
            DmapType::DMAP => {
                self.parse_record()?;
                DmapType::DMAP
            },
            _              => self.read_data(data_type)?
        };

        Ok(RawDmapScalar{data, name, mode})
    }

    fn read_data(&mut self, data_type: DmapType) -> Result<DmapType> {
        if self.cursor.position() > self.cursor.get_ref().len() as u64 {
            return Err(DmapError::new("READ DATA: Cursor extends out of buffer. Data is likely corrupted".to_string()))
        }
        if self.cursor.get_ref().len() as u64 - self.cursor.position() < data_type.get_num_bytes() {
            return Err(DmapError::new("READ DATA: Byte offsets into buffer are not properly aligned. \
            Data is likely corrupted".to_string()))
        }

        let position = self.cursor.position() as usize;
        let data_size = data_type.get_num_bytes() as usize;
        let data: &[u8] = &self.cursor.get_mut()[position..position+data_size];
        self.cursor.set_position({position + data_size} as u64);

        match data_type {
            DmapType::DMAP        => {
                match self.parse_record() {
                    Ok(_) => Ok(DmapType::DMAP),
                    Err(e) => Err(e)
                }
            },
            DmapType::UCHAR {..}  => {
                Ok(DmapType::UCHAR(<u8>::unpack_from_be(&mut data.clone())
                    .map_err(|_| DmapError::new("READ DATA: Unable to interpret data".to_string()))?))
            },
            DmapType::CHAR {..}   => {
                Ok(DmapType::CHAR(<i8>::unpack_from_be(&mut data.clone())
                    .map_err(|_| DmapError::new("READ DATA: Unable to interpret data".to_string()))?))
            },
            DmapType::SHORT {..}  => {
                Ok(DmapType::SHORT(<i16>::unpack_from_be(&mut data.clone())
                    .map_err(|_| DmapError::new("READ DATA: Unable to interpret data".to_string()))?))
            },
            DmapType::USHORT {..} => {
                Ok(DmapType::USHORT(<u16>::unpack_from_be(&mut data.clone())
                    .map_err(|_| DmapError::new("READ DATA: Unable to interpret data".to_string()))?))
            }
            DmapType::INT {..}    => {
                Ok(DmapType::INT(<i32>::unpack_from_be(&mut data.clone())
                    .map_err(|_| DmapError::new("READ DATA: Unable to interpret data".to_string()))?))
            },
            DmapType::UINT {..}   => {
                Ok(DmapType::UINT(<u32>::unpack_from_be(&mut data.clone())
                    .map_err(|_| DmapError::new("READ DATA: Unable to interpret data".to_string()))?))
            }
            DmapType::LONG {..}   => {
                Ok(DmapType::LONG(<i64>::unpack_from_be(&mut data.clone())
                    .map_err(|_| DmapError::new("READ DATA: Unable to interpret data".to_string()))?))
            }
            DmapType::ULONG {..}  => {
                Ok(DmapType::ULONG(<u64>::unpack_from_be(&mut data.clone())
                    .map_err(|_| DmapError::new("READ DATA: Unable to interpret data".to_string()))?))
            }
            DmapType::FLOAT {..}  => {
                Ok(DmapType::FLOAT(<f32>::unpack_from_be(&mut data.clone())
                    .map_err(|_| DmapError::new("READ DATA: Unable to interpret data".to_string()))?))
            }
            DmapType::DOUBLE {..} => {
                Ok(DmapType::DOUBLE(<f64>::unpack_from_be(&mut data.clone())
                    .map_err(|_| DmapError::new("READ DATA: Unable to interpret data".to_string()))?))
            }
            DmapType::STRING {..} => {
                let mut byte_counter = 0;
                while self.cursor.get_ref()[position + byte_counter] != 0 {
                    byte_counter += 1;
                    if position + byte_counter >= self.cursor.get_ref().len() {
                        return Err(DmapError::new("READ DATA: String is improperly terminated. \
                        Dmap record is corrupted".to_string()))
                    }
                }
                let mut data = String::from_utf8(self.cursor.get_ref().clone())
                    .map_err(|_| DmapError::new("READ DATA: Unable to interpret data".to_string()))?;
                self.cursor.set_position({position + byte_counter + 1} as u64);
                Ok(DmapType::STRING(data))
            }
        }
    }

}