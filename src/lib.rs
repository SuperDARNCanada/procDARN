use bytemuck;
use bytemuck::PodCastError;
use std::error::Error;
use std::fmt;
use std::fmt::Display;
use std::fs::File;
use std::io::{Cursor, Read, Write};
use std::path::Path;

type Result<T> = std::result::Result<T, DmapError>;

#[derive(Debug, Clone)]
pub enum DmapError {
    Parse(String, Vec<u8>),
    BadVal(String, DmapType),
    Message(String),
    CastError(String, PodCastError),
}

impl Error for DmapError {}

impl Display for DmapError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            DmapError::Message(msg) => write!(f, "{}", msg),
            DmapError::BadVal(msg, val) => write!(f, "{}: {:?}", msg, val),
            DmapError::Parse(msg, val) => write!(f, "{}: {:?}", msg, val),
            DmapError::CastError(msg, err) => write!(f, "{}: {}", msg, err.to_string()),
        }
    }
}

#[derive(Debug, Clone)]
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
    fn all_keys() -> Vec<i8> {
        vec![0, 1, 2, 3, 4, 8, 9, 10, 16, 17, 18, 19]
    }

    fn get_num_bytes(&self) -> u64 {
        match self {
            DmapType::CHAR { .. } => 1,
            DmapType::SHORT { .. } => 2,
            DmapType::INT { .. } => 4,
            DmapType::FLOAT { .. } => 4,
            DmapType::DOUBLE { .. } => 8,
            DmapType::LONG { .. } => 8,
            DmapType::UCHAR { .. } => 1,
            DmapType::USHORT { .. } => 2,
            DmapType::UINT { .. } => 4,
            DmapType::ULONG { .. } => 8,
            _ => 0,
        }
    }

    fn get_type_from_key(key: i8) -> Result<DmapType> {
        match key {
            0 => Ok(DmapType::DMAP),
            1 => Ok(DmapType::CHAR(0)),
            2 => Ok(DmapType::SHORT(0)),
            3 => Ok(DmapType::INT(0)),
            4 => Ok(DmapType::FLOAT(0.0)),
            8 => Ok(DmapType::DOUBLE(0.0)),
            9 => Ok(DmapType::STRING("".to_string())),
            10 => Ok(DmapType::LONG(0)),
            16 => Ok(DmapType::UCHAR(0)),
            17 => Ok(DmapType::USHORT(0)),
            18 => Ok(DmapType::UINT(0)),
            19 => Ok(DmapType::ULONG(0)),
            _ => Err(DmapError::Message(format!(
                "Invalid key for DMAP type: {}",
                key
            ))),
        }
    }

    fn to_bytes(&self) -> Vec<u8> {
        match self {
            DmapType::DMAP => vec![],
            DmapType::CHAR(x) => bytemuck::bytes_of(x).to_vec(),
            DmapType::UCHAR(x) => bytemuck::bytes_of(x).to_vec(),
            DmapType::SHORT(x) => bytemuck::bytes_of(x).to_vec(),
            DmapType::USHORT(x) => bytemuck::bytes_of(x).to_vec(),
            DmapType::INT(x) => bytemuck::bytes_of(x).to_vec(),
            DmapType::UINT(x) => bytemuck::bytes_of(x).to_vec(),
            DmapType::LONG(x) => bytemuck::bytes_of(x).to_vec(),
            DmapType::ULONG(x) => bytemuck::bytes_of(x).to_vec(),
            DmapType::FLOAT(x) => bytemuck::bytes_of(x).to_vec(),
            DmapType::DOUBLE(x) => bytemuck::bytes_of(x).to_vec(),
            DmapType::STRING(x) => {
                let mut bytes = vec![];
                bytes.append(&mut x.as_bytes().to_vec());
                bytes.push(0); // Rust String not null-terminated
                bytes
            }
        }
    }
}

impl Display for DmapType {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            DmapType::DMAP => write!(f, "DMAP"),
            DmapType::CHAR(x) => write!(f, "{}", x),
            DmapType::SHORT(x) => write!(f, "{}", x),
            DmapType::INT(x) => write!(f, "{}", x),
            DmapType::FLOAT(x) => write!(f, "{}", x),
            DmapType::DOUBLE(x) => write!(f, "{}", x),
            DmapType::STRING(x) => write!(f, "{:?}", x),
            DmapType::LONG(x) => write!(f, "{}", x),
            DmapType::UCHAR(x) => write!(f, "{}", x),
            DmapType::USHORT(x) => write!(f, "{}", x),
            DmapType::UINT(x) => write!(f, "{}", x),
            DmapType::ULONG(x) => write!(f, "{}", x),
        }
    }
}

#[derive(Debug)]
struct RawDmapScalar {
    data: DmapType,
    name: String,
    mode: i8,
}

impl RawDmapScalar {
    fn to_bytes(&self) -> Vec<u8> {
        let mut bytes: Vec<u8> = vec![];
        bytes.append(&mut DmapType::STRING(self.name.clone()).to_bytes());
        bytes.append(&mut DmapType::CHAR(self.mode).to_bytes());
        bytes.append(&mut self.data.to_bytes());
        bytes
    }
}

#[derive(Debug)]
struct RawDmapArray {
    name: String,
    mode: i8,
    _arr_dimensions: Vec<i32>,
    data: Vec<DmapType>,
}

impl RawDmapArray {
    fn to_bytes(&self) -> Vec<u8> {
        let mut bytes: Vec<u8> = vec![];
        bytes.append(&mut DmapType::STRING(self.name.clone()).to_bytes());
        bytes.append(&mut DmapType::CHAR(self.mode).to_bytes());
        for val in self.data.clone() {
            bytes.append(&mut val.to_bytes());
        }
        bytes
    }
}

#[derive(Debug)]
pub struct RawDmapRecord {
    num_scalars: i32,
    num_arrays: i32,
    scalars: Vec<RawDmapScalar>,
    arrays: Vec<RawDmapArray>,
}

impl RawDmapRecord {
    fn to_bytes(&self) -> Vec<u8> {
        let mut container: Vec<u8> = vec![];
        let code = 65537; // No idea why this is what it is, copied from backscatter

        let mut data_bytes: Vec<u8> = vec![];
        for scalar in &self.scalars {
            data_bytes.append(&mut scalar.to_bytes());
        }
        for array in &self.arrays {
            data_bytes.append(&mut array.to_bytes());
        }

        container.append(&mut DmapType::INT(code).to_bytes());
        container.append(&mut DmapType::INT(data_bytes.len() as i32).to_bytes());
        container.append(&mut DmapType::INT(self.num_scalars).to_bytes());
        container.append(&mut DmapType::INT(self.num_arrays).to_bytes());
        container
    }
}

#[derive(Debug)]
pub struct RawDmapRead {
    cursor: Cursor<Vec<u8>>,
    pub dmap_records: Vec<RawDmapRecord>,
}

impl RawDmapRead {
    pub fn new(mut dmap_data: impl Read) -> Result<RawDmapRead> {
        let mut buffer: Vec<u8> = vec![];

        dmap_data
            .read_to_end(&mut buffer)
            .map_err(|_| DmapError::Message("Could not read data".to_string()))?;

        let cursor = Cursor::new(buffer);
        let mut dmap_read = RawDmapRead {
            cursor,
            dmap_records: vec![],
        };

        while dmap_read.cursor.position() < dmap_read.cursor.get_ref().len() as u64 {
            dmap_read.parse_record()?;
        }
        Ok(dmap_read)
    }

    fn parse_record(&mut self) -> Result<()> {
        let bytes_already_read = self.cursor.position();
        let _code = match self.read_data(DmapType::INT(0))? {
            DmapType::INT(i) => Ok(i),
            _ => Err(DmapError::Message("PARSE RECORD: Invalid code".to_string())),
        }?;
        let size = match self.read_data(DmapType::INT(0))? {
            DmapType::INT(i) => Ok(i),
            _ => Err(DmapError::Message("PARSE RECORD: Invalid size".to_string())),
        }?;

        // adding 8 bytes because code and size are part of the record.
        if size as u64
            > self.cursor.get_ref().len() as u64 - self.cursor.position()
                + 2 * DmapType::INT(0).get_num_bytes()
        {
            return Err(DmapError::Message(
                "PARSE RECORD: Integrity check shows record size bigger than \
                remaining buffer. Data is likely corrupted"
                    .to_string(),
            ));
        } else if size <= 0 {
            return Err(DmapError::Message(
                "PARSE RECORD: Integrity check shows record size <= 0. \
                Data is likely corrupted"
                    .to_string(),
            ));
        }

        let num_scalars = match self.read_data(DmapType::INT(0))? {
            DmapType::INT(i) => Ok(i),
            _ => Err(DmapError::Message(
                "PARSE RECORD: Invalid number of scalars".to_string(),
            )),
        }?;
        let num_arrays = match self.read_data(DmapType::INT(0))? {
            DmapType::INT(i) => Ok(i),
            _ => Err(DmapError::Message(
                "PARSE RECORD: Invalid number of arrays".to_string(),
            )),
        }?;
        if num_scalars <= 0 {
            return Err(DmapError::Message(
                "PARSE RECORD: Number of scalers is 0 or negative.".to_string(),
            ));
        } else if num_arrays <= 0 {
            return Err(DmapError::Message(
                "PARSE RECORD: Number of arrays is 0 or negative.".to_string(),
            ));
        } else if num_scalars + num_arrays > size {
            return Err(DmapError::Message(
                "PARSE RECORD: Invalid number of record elements. \
                Array or scaler field is likely corrupted."
                    .to_string(),
            ));
        }

        let mut scalars: Vec<RawDmapScalar> = vec![];
        for _ in 0..num_scalars {
            scalars.push(self.parse_scalar()?);
        }

        let mut arrays: Vec<RawDmapArray> = vec![];
        for _ in 0..num_arrays {
            arrays.push(self.parse_array(size)?);
        }

        if self.cursor.position() - bytes_already_read != size as u64 {
            return Err(DmapError::Message(format!(
                "PARSE RECORD: Bytes read {} does not match the records size field {}",
                self.cursor.position() - bytes_already_read,
                size
            )));
        }

        self.dmap_records.push(RawDmapRecord {
            num_scalars,
            scalars,
            num_arrays,
            arrays,
        });
        Ok(())
    }

    fn parse_scalar(&mut self) -> Result<RawDmapScalar> {
        let mode = 6;
        let name = match self.read_data(DmapType::STRING("".to_string()))? {
            DmapType::STRING(s) => Ok(s),
            _ => Err(DmapError::Message(
                "PARSE SCALAR: Invalid scalar name".to_string(),
            )),
        }?;
        let data_type_key = match self.read_data(DmapType::CHAR(0))? {
            DmapType::CHAR(c) => Ok(c),
            _ => Err(DmapError::Message(
                "PARSE SCALAR: Invalid data type".to_string(),
            )),
        }?;

        if !DmapType::all_keys().contains(&data_type_key) {
            return Err(DmapError::BadVal(
                "PARSE SCALAR: Data type is corrupted. Record is likely \
                corrupted"
                    .to_string(),
                DmapType::CHAR(data_type_key),
            ));
        }

        let data_type = DmapType::get_type_from_key(data_type_key)?;

        let data = match data_type {
            DmapType::DMAP => {
                self.parse_record()?;
                DmapType::DMAP
            }
            _ => self.read_data(data_type)?,
        };

        Ok(RawDmapScalar { data, name, mode })
    }

    fn parse_array(&mut self, record_size: i32) -> Result<RawDmapArray> {
        let mode = 7;
        let name = match self.read_data(DmapType::STRING("".to_string()))? {
            DmapType::STRING(s) => Ok(s),
            _ => Err(DmapError::Message(
                "PARSE ARRAY: Invalid array name".to_string(),
            )),
        }?;
        let data_type_key = match self.read_data(DmapType::CHAR(0))? {
            DmapType::CHAR(c) => Ok(c),
            _ => Err(DmapError::Message(
                "PARSE ARRAY: Invalid data type".to_string(),
            )),
        }?;

        if !DmapType::all_keys().contains(&data_type_key) {
            return Err(DmapError::Message(
                "PARSE ARRAY: Data type is corrupted. Record is likely \
                corrupted"
                    .to_string(),
            ));
        }

        let data_type = DmapType::get_type_from_key(data_type_key)?;

        let array_dimension = match self.read_data(DmapType::INT(0))? {
            DmapType::INT(i) => Ok(i),
            _ => Err(DmapError::Message(
                "PARSE ARRAY: Invalid array dimension".to_string(),
            )),
        }?;

        if array_dimension > record_size {
            return Err(DmapError::Message(
                "PARSE ARRAY: Parsed # of array dimensions are larger \
                than record size. Record is likely corrupted"
                    .to_string(),
            ));
        } else if array_dimension <= 0 {
            return Err(DmapError::Message(
                "PARSE ARRAY: Parsed # of array dimensions are zero or \
                negative. Record is likely corrupted"
                    .to_string(),
            ));
        }

        let mut dimensions: Vec<i32> = vec![];
        let mut total_elements = 1;
        for _ in 0..array_dimension {
            let dim = match self.read_data(DmapType::INT(0))? {
                DmapType::INT(val) => Ok(val),
                _ => Err(DmapError::Message(
                    "PARSE ARRAY: Array dimensions could not be parsed".to_string(),
                )),
            }?;
            if dim <= 0 {
                return Err(DmapError::Message(
                    "PARSE ARRAY: Array dimension is zero or negative. \
                    Record is likely corrupted"
                        .to_string(),
                ));
            } else if dim > record_size {
                return Err(DmapError::Message(
                    "PARSE ARRAY: Array dimension exceeds record size".to_string(),
                ));
            }
            dimensions.push(dim);
            total_elements = total_elements * dim;
        }

        if total_elements > record_size {
            return Err(DmapError::Message(
                "PARSE ARRAY: Total array elements > record size.".to_string(),
            ));
        } else if total_elements * data_type.get_num_bytes() as i32 > record_size {
            return Err(DmapError::Message(
                "PARSE ARRAY: Array size exceeds record size. Data is \
                likely corrupted"
                    .to_string(),
            ));
        }
        let mut data = vec![];
        for _ in 0..total_elements {
            data.push(self.read_data(data_type.clone())?);
        }
        Ok(RawDmapArray {
            name,
            mode,
            _arr_dimensions: dimensions,
            data,
        })
    }

    fn read_data(&mut self, data_type: DmapType) -> Result<DmapType> {
        if self.cursor.position() > self.cursor.get_ref().len() as u64 {
            return Err(DmapError::Message(
                "READ DATA: Cursor extends out of buffer. Data is likely corrupted".to_string(),
            ));
        }
        if self.cursor.get_ref().len() as u64 - self.cursor.position() < data_type.get_num_bytes() {
            return Err(DmapError::Message(
                "READ DATA: Byte offsets into buffer are not properly aligned. \
            Data is likely corrupted"
                    .to_string(),
            ));
        }

        let position = self.cursor.position() as usize;
        let mut data_size = data_type.get_num_bytes() as usize;
        let data: &[u8] = &self.cursor.get_mut()[position..position + data_size];
        let parsed_data = match data_type {
            DmapType::DMAP => self.parse_record().map(|_| DmapType::DMAP)?,
            DmapType::UCHAR { .. } => DmapType::UCHAR(data[0]),
            DmapType::CHAR { .. } => {
                DmapType::CHAR(*bytemuck::try_from_bytes::<i8>(data).map_err(|_| {
                    DmapError::Message("READ DATA: Unable to interpret char".to_string())
                })?)
            }
            DmapType::SHORT { .. } => {
                DmapType::SHORT(bytemuck::try_pod_read_unaligned::<i16>(data).map_err(|e| {
                    DmapError::CastError("READ DATA: Unable to interpret short".to_string(), e)
                })?)
            }
            DmapType::USHORT { .. } => {
                DmapType::USHORT(*bytemuck::try_from_bytes::<u16>(data).map_err(|e| {
                    DmapError::CastError("READ DATA: Unable to interpret ushort".to_string(), e)
                })?)
            }
            DmapType::INT { .. } => {
                DmapType::INT(bytemuck::try_pod_read_unaligned::<i32>(data).map_err(|e| {
                    DmapError::CastError("READ DATA: Unable to interpret int".to_string(), e)
                })?)
            }
            DmapType::UINT { .. } => {
                DmapType::UINT(*bytemuck::try_from_bytes::<u32>(data).map_err(|_| {
                    DmapError::Message("READ DATA: Unable to interpret uint".to_string())
                })?)
            }
            DmapType::LONG { .. } => {
                DmapType::LONG(*bytemuck::try_from_bytes::<i64>(data).map_err(|_| {
                    DmapError::Message("READ DATA: Unable to interpret long".to_string())
                })?)
            }
            DmapType::ULONG { .. } => {
                DmapType::ULONG(*bytemuck::try_from_bytes::<u64>(data).map_err(|_| {
                    DmapError::Message("READ DATA: Unable to interpret ulong".to_string())
                })?)
            }
            DmapType::FLOAT { .. } => {
                DmapType::FLOAT(bytemuck::try_pod_read_unaligned::<f32>(data).map_err(|_| {
                    DmapError::Message("READ DATA: Unable to interpret float".to_string())
                })?)
            }
            DmapType::DOUBLE { .. } => {
                DmapType::DOUBLE(bytemuck::try_pod_read_unaligned::<f64>(data).map_err(|_| {
                    DmapError::Message("READ DATA: Unable to interpret double".to_string())
                })?)
            }
            DmapType::STRING { .. } => {
                let mut byte_counter = 0;
                let stream = self.cursor.get_ref();
                while stream[position + byte_counter] != 0 {
                    byte_counter += 1;
                    if position + byte_counter >= stream.len() {
                        return Err(DmapError::Message(
                            "READ DATA: String is improperly terminated. \
                        Dmap record is corrupted"
                                .to_string(),
                        ));
                    }
                }
                let data = String::from_utf8(stream[position..position + byte_counter].to_owned())
                    .map_err(|_| {
                        DmapError::Message("READ DATA: Unable to interpret string".to_string())
                    })?;
                data_size = byte_counter + 1;
                DmapType::STRING(data)
            }
        };
        self.cursor.set_position({ position + data_size } as u64);

        Ok(parsed_data)
    }
}

pub fn to_file<P: AsRef<Path>>(path: P, dmap_records: Vec<RawDmapRecord>) -> std::io::Result<()> {
    // let writer = RawDmapWrite{dmap_records, stream: vec![]};
    let mut stream = vec![];
    for rec in dmap_records {
        stream.append(&mut rec.to_bytes());
    }
    let mut file = File::create(path)?;
    file.write_all(&stream)?;
    Ok(())
}
