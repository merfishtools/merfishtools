use custom_derive::custom_derive;
use newtype_derive::*;

custom_derive! {
    #[derive(Copy, Clone, PartialEq, PartialOrd, Eq, Ord, Hash,
        NewtypeFrom, NewtypeAdd, NewtypeBinary,
        NewtypeDeref, NewtypeDerefMut)]
    pub struct Barcode(pub u16);
}

impl std::fmt::Debug for Barcode {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> Result<(), std::fmt::Error> {
        f.write_fmt(format_args!("{:016b}", self))
    }
}