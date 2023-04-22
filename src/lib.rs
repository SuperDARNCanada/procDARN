use std::io::{Read};
use ndarray::prelude::*;
use structure;

#[derive(Debug)]
struct DmapType {
    fmt: char
}

impl DmapType {
    fn get_num_bytes(self) -> usize {
        match self.fmt {
            'c' => 1,
            'B' => 1,
            'h' => 2,
            'H' => 2,
            'i' => 4,
            'I' => 4,
            'q' => 8,
            'Q' => 8,
            'f' => 4,
            'd' => 8,
            _ => 0
        }
    }
}

#[derive(Debug)]
pub enum DmapTypes {
    DMAP = 0,
    CHAR = 1,
    SHORT = 2,
    INT = 3,
    FLOAT = 4,
    DOUBLE = 8,
    STRING = 9,
    LONG = 10,
    UCHAR = 16,
    USHORT = 17,
    UINT = 18,
    ULONG = 19
}

impl DmapTypes {
    fn all_keys() -> Vec<usize> {
        vec![0, 1, 2, 3, 4, 8, 9, 10, 16, 17, 18, 19]
    }
}

struct RawDmapScalar {
    dmap_type: DmapTypes,
    name: str,
    mode: str,
    data: str,
    data_type_fmt: char
}

struct RawDmapArray {
    dmap_type: DmapTypes,
    name: str,
    mode: str,
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
    dmap_records: Vec<RawDmapRecord>
}


impl RawDmapRead {

    fn new(dmap_data: &mut &impl Read) -> RawDmapRead {
        let mut cursor: usize = 0;
        let mut buffer: Vec<u8>;
        dmap_data.read_to_end(&mut buffer)?;

        // TODO: Test initial data integrity

        while cursor < buffer.len() {
            parse_record(cursor, buffer);
        }
        Ok(())
    }

    fn read_data(&self, buffer: Vec<u8>, data_type_fmt: char) {
        if self.cursor > buffer.len() {
            println!("READ DATA: Cursor extends out of buffer. Data is likely corrupted")
            // TODO: Raise DmapDataError
        }
        if buffer.len() - self.cursor < self.get_num_bytes(data_type_fmt) {
            println!("READ DATA: Byte offsets into buffer are not properly aligned. \
            Data is likely corrupted")
            // TODO: Raise DmapDataError
        }

        if data_type_fmt == DmapTypes.DMAP {
            self.parse_record()
        } else if data_type_fmt == 'c' {
            let mut data = buffer[self.cursor];
            self.cursor += get_num_bytes(data_type_fmt);
        } else if data_type_fmt != 's' {
            let s = structure!(data_type_fmt);
            let mut data = s.unpack_from(&mut buffer, self.cursor)?;
            self.cursor += get_num_bytes(data_type_fmt);
        } else {
            let mut byte_counter = 0;
            while buffer[self.cursor + byte_counter] != 0 {
                byte_counter += 1;
                if self.cursor + byte_counter >= buffer.len() {
                    println!("READ DATA: String is improperly terminated. Dmap record is corrupted")
                    // TODO: Raise DmapDataError
                }
            }
            let char_count = String::from("{}s", byte_counter);
            let s = structure!(char_count);
            let mut data = s.unpack_from(&mut buffer, self.cursor)?;
            data = (data[0],);
            self.cursor += byte_counter + 1;
        }

        if data_type_fmt == 'c' { data }
        else { data[0] }
    }

}