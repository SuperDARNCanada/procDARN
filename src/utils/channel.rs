use crate::error::ProcdarnError;

pub fn set_stereo_channel(channel_char: char) -> Result<i32, ProcdarnError> {
    match channel_char {
        'a' => Ok(1),
        'b' => Ok(2),
        _ => Err(ProcdarnError::Channel(format!(
            "Invalid stereo channel {}",
            channel_char
        ))),
    }
}

pub fn set_fix_channel(channel_char: char) -> Result<i32, ProcdarnError> {
    match channel_char {
        'a' => Ok(1),
        'b' => Ok(2),
        'c' => Ok(3),
        'd' => Ok(4),
        _ => Err(ProcdarnError::Channel(format!(
            "Invalid fix channel {}",
            channel_char
        ))),
    }
}
