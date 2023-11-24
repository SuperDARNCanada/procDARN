use crate::error::BackscatterError;

pub fn set_stereo_channel(channel_char: char) -> Result<i32, BackscatterError> {
    match channel_char {
        'a' => Ok(1),
        'b' => Ok(2),
        _ => Err(BackscatterError::new(
            format!("Invalid stereo channel {}", channel_char).as_str(),
        )),
    }
}

pub fn set_fix_channel(channel_char: char) -> Result<i32, BackscatterError> {
    match channel_char {
        'a' => Ok(1),
        'b' => Ok(2),
        'c' => Ok(3),
        'd' => Ok(4),
        _ => Err(BackscatterError::new(
            format!("Invalid fix channel {}", channel_char).as_str(),
        )),
    }
}
