use std::{fmt, mem};
use std::any::TypeId;
use std::intrinsics::{size_of, size_of_val};
use std::io::{Read};
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
    CHAR,
    SHORT,
    INT,
    FLOAT,
    DOUBLE,
    STRING,
    LONG,
    UCHAR,
    USHORT,
    UINT,
    ULONG,
}

impl DmapType {
    fn all_keys() -> Vec<usize> { vec![0, 1, 2, 3, 4, 8, 9, 10, 16, 17, 18, 19] }

    fn get_num_bytes(self) -> usize {
        match self {
            DmapType::CHAR   => 1,
            DmapType::SHORT  => 2,
            DmapType::INT    => 4,
            DmapType::FLOAT  => 4,
            DmapType::DOUBLE => 8,
            DmapType::LONG   => 8,
            DmapType::UCHAR  => 1,
            DmapType::USHORT => 2,
            DmapType::UINT   => 4,
            DmapType::ULONG  => 8,
            _                => 0
        }
    }

    fn get_fmt(self) -> TypeId {
        match self {
            DmapType::CHAR   => i8,
            DmapType::SHORT  => i16,
            DmapType::INT    => i32,
            DmapType::FLOAT  => f32,
            DmapType::DOUBLE => f64,
            DmapType::STRING => u8,
            DmapType::LONG   => i64,
            DmapType::UCHAR  => u8,
            DmapType::USHORT => u16,
            DmapType::UINT   => u32,
            DmapType::ULONG  => u64,
            DmapType::DMAP   => u8,
        }
    }

    fn get_type_from_fmt(fmt: char) -> Result<DmapType> {
        match fmt {
            'c' => Ok(DmapType::CHAR),
            'h' => Ok(DmapType::SHORT),
            'i' => Ok(DmapType::INT),
            'f' => Ok(DmapType::FLOAT),
            'd' => Ok(DmapType::DOUBLE),
            's' => Ok(DmapType::STRING),
            'q' => Ok(DmapType::LONG),
            'B' => Ok(DmapType::UCHAR),
            'H' => Ok(DmapType::USHORT),
            'I' => Ok(DmapType::UINT),
            'Q' => Ok(DmapType::ULONG),
            _   => Err(DmapError::new(format!("Invalid format for DMAP type: {}", fmt)))
        }
    }

    fn get_type_from_key(key: usize) -> Result<DmapType> {
        match key {
            0  => Ok(DmapType::DMAP),
            1  => Ok(DmapType::CHAR),
            2  => Ok(DmapType::SHORT),
            3  => Ok(DmapType::INT),
            4  => Ok(DmapType::FLOAT),
            8  => Ok(DmapType::DOUBLE),
            9  => Ok(DmapType::STRING),
            10 => Ok(DmapType::LONG),
            16 => Ok(DmapType::UCHAR),
            17 => Ok(DmapType::USHORT),
            18 => Ok(DmapType::UINT),
            19 => Ok(DmapType::ULONG),
            _  => Err(DmapError::new(format!("Invalid key for DMAP type: {}", key)))
        }
    }
}

struct RawDmapScalar {
    dmap_type: DmapType,
    name: Box<str>,
    mode: Box<str>,
    data: Box<str>,
    data_type_fmt: char
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
    num_scalars: u32,
    num_arrays: u32,
    scalars: Vec<RawDmapScalar>,
    arrays: Vec<RawDmapArray>
}

struct RawDmapRead {
    cursor: usize,
    dmap_records: Vec<RawDmapRecord>,
    buffer: Vec<u8>
}

impl RawDmapRead {

    fn new(dmap_data: &mut &impl Read) -> Result<RawDmapRead> {
        let mut buffer: Vec<u8>;

        dmap_data.read_to_end(&mut buffer)
            .ok_or(DmapError::new("Could not read data".to_string()));

        let mut dmap_read = RawDmapRead{cursor: 0, dmap_records: vec![], buffer};

        // TODO: Test initial data integrity

        while dmap_read.cursor < dmap_read.buffer.len() {
            dmap_read.parse_record();
        }
        Ok(dmap_read)
    }

    fn parse_record(&self) -> Result<()> {
        let bytes_already_read = self.cursor;
        let code = self.read_data(DmapType::INT)?;
        let size = self.read_data(DmapType::INT)?;

        // adding 8 bytes because code and size are part of the record.
        if size as usize > self.buffer.len() - self.cursor + 2 * DmapType::INT.get_num_bytes() {
            return Err(DmapError::new("PARSE RECORD: Integrity check shows record size bigger than \
                remaining buffer. Data is likely corrupted".to_string()))
        }
        else if size <= 0 {
            return Err(DmapError::new("PARSE RECORD: Integrity check shows record size <= 0. \
                Data is likely corrupted".to_string()))
        }

        let num_scalars = self.read_data(DmapType::INT)?;
        let num_arrays = self.read_data(DmapType::INT)?;

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

        let scalars: Vec<RawDmapScalar>;
        for n in 0..num_scalars {
            scalars.push(self.parse_scalar());
        }
        let arrays: Vec<RawDmapArray>;
        for i in 0..num_arrays {
            arrays.push(self.parse_array(size));
        }

        if self.cursor - bytes_already_read != size as usize {
            return Err(DmapError::new(format!(
                "PARSE RECORD: Bytes read {} does not match the records size field {}",
                self.cursor - bytes_already_read, size)))
        }

        self.dmap_records.push(RawDmapRecord{num_scalars: num_scalars.into(), scalars,
            num_arrays: num_arrays.into(), arrays});
        Ok(())
    }

    fn read_data(&self, data_type: DmapType) -> Result<u8> {
        if self.cursor > self.buffer.len() {
            println!("READ DATA: Cursor extends out of buffer. Data is likely corrupted")
            // TODO: Raise DmapDataError
        }
        if self.buffer.len() - self.cursor < data_type.get_num_bytes() {
            println!("READ DATA: Byte offsets into buffer are not properly aligned. \
            Data is likely corrupted")
            // TODO: Raise DmapDataError
        }

        match data_type {
            DmapType::DMAP => {
                self.parse_record()?;
                Ok(())
            }
            DmapType::CHAR => {
                let mut data = self.buffer[self.cursor];
                self.cursor += data_type.get_num_bytes();
                data
            }
            DmapType::STRING => {
                let mut byte_counter = 0;
                while self.buffer[self.cursor + byte_counter] != 0 {
                    byte_counter += 1;
                    if self.cursor + byte_counter >= self.buffer.len() {
                        println!("READ DATA: String is improperly terminated. \
                        Dmap record is corrupted")
                        // TODO: Raise DmapDataError
                    }
                }
                // let char_count = format!("{}s", byte_counter);
                // let s = structure::structure!(char_count.as_str());
                // let mut data = s.unpack_from(&mut self.buffer, self.cursor)?;
                let mut data = <>::unpack_from_be(&mut self.buffer[])?;
                data = (data[0],);
                self.cursor += byte_counter + 1;
                data[0]     // structure.unpack_from returns a tuple. [0] is the actual data
            }
            _ => {
                let mut data = <>::unpack_from_be(&mut self.buffer, self.cursor)?;
                self.cursor += data_type.get_num_bytes();
                data[0]     // structure.unpack_from returns a tuple. [0] is the actual data
            }
        }
    }

}