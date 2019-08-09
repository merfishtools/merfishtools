use std::fs;
use std::io;
use std::io::Read;
use std::path::Path;

use bincode;
use bit_vec::BitVec;
use byteorder::{ByteOrder, NativeEndian};
use failure::Error;

use crate::io::counts::{CommonRecord, FromRecord};
use crate::io::merfishdata::{MerfishRecord, Readout};

/// Header of a binary merfish file.
#[derive(Serialize, Deserialize, Debug)]
pub struct Header {
    version: u8,
    #[serde(rename = "isCorrupt")]
    is_corrupt: u8,
    #[serde(rename = "numEntries")]
    num_entries: u32,
    #[serde(rename = "headerLength")]
    header_length: u32,
}

impl Header {
    pub fn version(&self) -> u8 {
        self.version
    }
    pub fn is_corrupt(&self) -> u8 {
        self.is_corrupt
    }
    pub fn num_entries(&self) -> u32 {
        self.num_entries
    }
    /// Length of the line giving the field names of a record
    pub fn header_length(&self) -> u32 {
        self.header_length
    }
}

/// A `Record` represents a single entry of a merfish binary file.
#[derive(Serialize, Deserialize, Debug, PartialEq, Clone)]
pub struct BinaryRecord {
    /// An integer representation of the 16-bit barcode associated with each RNA.
    pub barcode: u64,
    /// The specific codebook entry to which this barcode corresponds, e.g. '1' corresponds to the first RNA listed in the codebook ('Blank-1'), '2' corresponds to the second ('Blank-10'), etc.
    pub barcode_id: u16,
    /// The id associated with the field-of-view in which this RNA was imaged.
    pub fov_id: u16,
    /// The sum of the normalized intensity associated with each pixel assigned to a given RNA.
    pub total_magnitude: f32,
    /// The center of the pixels associated with a given RNA in the coordinate system of each field-of-view in units of pixels.
    pub pixel_centroid: [u16; 2],
    /// The center of an RNA, weighted by the intensity of each pixel, in the coordinate system of each field-of-view in units of pixels.
    pub weighted_pixel_centroid: [f32; 2],
    /// The center of an RNA in the absolute coordinate system of the sample in microns. This quantity is calculated from the weighted_pixel_centroid.
    pub abs_position: [f32; 2],
    /// The number of pixels associated with an RNA.
    pub area: u16,
    /// The normalized intensity of each pixel for each bit in the barcode, average across all pixels assigned to an RNA.
    pub pixel_trace_mean: [f32; 16],
    /// The standard deviation of the normalized intensity for each bit in the barcode taken across all pixels assigned to an RNA.
    pub pixel_trace_std: [f32; 16],
    /// A boolean flag representing whether or not (TRUE/FALSE) the measured pixel trace for an RNA matched exactly to a barcode or required error correction at one bit.
    pub is_exact: u8,
    /// The bit at which error correction was applied. (0 if it was not applied).
    pub error_bit: u8,
    /// A boolean representing the type of error (0 corresponds to a '0' to '1' error; 1 corresponds to a '1' to '0' error).
    pub error_dir: u8,
    /// The average Euclidean distance between the pixel trace for a given RNA and the barcode to which it was matched.
    pub av_distance: f32,
    /// An integer representing the cell to which an RNA was assigned. This value corresponds to the index in the corresponding cellBoundaries.mat file.
    #[serde(rename = "cellID")]
    pub cell_id: u32,
    /// A boolean flag representing whether or not (TRUE/FALSE) an RNA was found within the nucleus of the cell.
    #[serde(rename = "inNucleus")]
    pub in_nucleus: u8,
    /// The distance (in microns) from the RNA to the closest edge of the nucleus.
    #[serde(rename = "distNucleus")]
    pub dist_nucleus: f64,
    /// The distance (in microns) from the RNA to the closest boundary of the cell.
    #[serde(rename = "distPeriphery")]
    pub dist_periphery: f64,
}

impl MerfishRecord for BinaryRecord {
    fn cell_id(&self) -> u32 {
        self.cell_id
    }

    fn cell_name(&self) -> String {
        self.cell_id.to_string()
    }

    fn cell_pos(&self) -> (f32, f32) {
        (self.abs_position[0], self.abs_position[1])
    }

    fn feature_id(&self) -> u16 {
        self.barcode_id
    }

    fn feature_name(&self) -> String {
        self.barcode_id.to_string()
    }

    fn hamming_dist(&self) -> u8 {
        assert!(
            self.is_exact <= 1,
            "unexpected value in field is_exact: {}",
            self.is_exact
        );
        1 - self.is_exact
    }

    fn error_mask(&self) -> u16 {
        if self.error_bit > 0 {
            1 << (self.error_bit as u16 - 1)
        } else {
            0b0000_0000_0000_0000
        }
    }

    fn codeword(&self) -> Option<u16> {
        Some(self.barcode as u16 ^ self.error_mask())
    }

    fn readout(&self) -> u16 {
        self.barcode as u16
    }

    fn readout_bitvec(&self) -> Readout {
        let mut buf = [0; 8];
        NativeEndian::write_u64(&mut buf, self.barcode);
        let mut readout = BitVec::with_capacity(16);
        (0..16).for_each(|i| readout.push(((self.barcode >> i) & 1) == 1));
        readout.truncate(16);
        if self.error_bit > 0 {
            // error_bit == 0 <=> there is no error.
            // i.e. the actually erroneous bit is `error_bit - 1`
            let bit = self.error_bit as usize - 1;
            let value = !readout.get(bit).unwrap();
            readout.set(bit, value);
        }
        readout
    }

    fn count(&self) -> usize {
        1
    }

    fn is_exact(&self) -> bool { self.is_exact != 0 }
}

#[derive(Debug, Fail)]
pub enum ReaderError {
    #[fail(display = "unsupported version: {}", version)]
    UnsupportedVersion { version: u8 },
    #[fail(
    display = "header is corrupt, i.e. might not have been written properly (while in append mode)"
    )]
    Corrupt,
}

pub struct BinaryReader<R: io::Read> {
    reader: io::BufReader<R>,
    header: Header,
}

impl<R: io::Read> BinaryReader<R> {
    /// Constructs a new `merfish::Reader<R>` where `R: io::Read`.
    ///
    /// This will read the header in order to check if the version matches.
    pub fn new(reader: R) -> Result<Self, Error> {
        let mut reader = io::BufReader::new(reader);
        let header: Result<Header, _> = bincode::deserialize_from(&mut reader);
        match header {
            Ok(ref header) if header.version != 1 => Err(ReaderError::UnsupportedVersion {
                version: header.version,
            })?,
            Ok(ref header) if header.is_corrupt != 0 => Err(ReaderError::Corrupt)?,
            Ok(header) => {
                // Read (discard) table header
                let mut header_buf = vec![0u8; header.header_length as usize];
                reader.read_exact(&mut header_buf)?;
                Ok(BinaryReader { reader, header })
            }
            Err(e) => Err(e)?,
        }
    }

    pub fn header(&self) -> &Header {
        &self.header
    }

    pub fn read(&mut self) -> Result<BinaryRecord, bincode::Error> {
        bincode::deserialize_from::<_, BinaryRecord>(&mut self.reader)
    }
}

impl<'a, R: io::Read + 'a> super::Reader<'a> for BinaryReader<R> {
    type Record = BinaryRecord;
    type Error = bincode::Error;
    type Iterator = BinaryRecordIterator<'a, R>;

    fn records(&'a mut self) -> Self::Iterator {
        BinaryRecordIterator { reader: self, i: 0 }
    }
}

impl BinaryReader<fs::File> {
    /// Constructs a new `merfish::Reader<fs::File>` from given file.
    ///
    /// This delegates work to `merfish::Reader::new`.
    pub fn from_file<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        fs::File::open(path).map(|r| BinaryReader::new(r).unwrap())
    }
}

pub struct BinaryRecordIterator<'a, R: io::Read + 'a> {
    reader: &'a mut BinaryReader<R>,
    i: u32,
}

impl<'a, R: io::Read> Iterator for BinaryRecordIterator<'a, R> {
    type Item = Result<BinaryRecord, bincode::Error>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.i >= self.reader.header.num_entries {
            None
        } else {
            self.i += 1;
            Some(self.reader.read())
        }
    }
}


#[cfg(test)]
mod tests {
    use std::io;

    use super::super::*;

    #[test]
    fn test_binary_records() {
        let data = b"\x01\x00\x01\x00\x00\x00\xad\x01\x00\x00\
        barcode,1  1,uint64,barcode_id,1  1,uint16,fov_id,1  1,uint16,\
        total_magnitude,1  1,single,pixel_centroid,1  2,uint16,weighted_pixel_centroid,1  2,\
        single,abs_position,1  2,single,area,1  1,uint16,pixel_trace_mean,1  16,\
        single,pixel_trace_std,1  16,single,is_exact,1  1,uint8,error_bit,1  1,uint8,\
        error_dir,1  1,uint8,av_distance,1  1,single,cellID,1  1,uint32,inNucleus,1  1,uint8,\
        distNucleus,1  1,double,distPeriphery,1  1,double\
        \x1c\x10\x00\x00\x00\x00\x00\x00\x01\x00\x00\x00\xdb\x1e\xb4A\xae\x00\xe9\x05\x14\x87.C\xe2\x17\xbdD\x9b\xe7\xaf\xc4\x10*8\xc5\x03\x00\xdf>\x1b;\x943\x06><\xee\xeb>\x8b\x97\x95>z\x05C?j\x9c\xe6;lw]=\x1fWp;\xa5\xbfz8\x00\x00\x00\x00 _@<\x87\xff\x8f<\x0c\x11\x88>3\xb1\xdc=_\x03\x8c<\xcd\xee\x99=\xe8\x8c[;AN2=v\x0b\x1d=#S\xd5<\xed\x8a<;\xcf\x81<;wYL=\xf3\xb8;;YN\xb18\x00\x00\x00\x00\x89\x03\xa1;Xh\xc7<zwD=u\xfa:=L\xe9O<\xb7\xe3\x1a=\x01\x00\x00\xcc&\xef>\x02\x00\x00\x00\x01\x00\x00\x00\x80a\xaf\x08@\x00\x00\x00\x809\xfc(@";

        let mut reader = binary::BinaryReader::new(io::Cursor::new(&data[..])).unwrap();
        let expected_record = binary::BinaryRecord {
            barcode: 4124,
            barcode_id: 1,
            fov_id: 0,
            total_magnitude: 22.515066,
            pixel_centroid: [174, 1513],
            weighted_pixel_centroid: [174.52765, 1512.7463],
            abs_position: [-1407.2377, -2946.629],
            area: 3,
            pixel_trace_mean: [
                0.0023688597,
                0.13105613,
                0.46080196,
                0.2921718,
                0.7618023,
                0.007037689,
                0.054068968,
                0.0036673022,
                0.000059783128,
                0.0,
                0.01174143,
                0.0175779,
                0.26575506,
                0.107759856,
                0.017091451,
                0.07516251,
            ],
            pixel_trace_std: [
                0.0033500735,
                0.04353166,
                0.038341008,
                0.02604062,
                0.002876933,
                0.0028763895,
                0.049890008,
                0.0028644174,
                0.00008454611,
                0.0,
                0.0049137515,
                0.024341747,
                0.047965504,
                0.04564901,
                0.0126899,
                0.037814822,
            ],
            is_exact: 1,
            error_bit: 0,
            error_dir: 0,
            av_distance: 0.46709287,
            cell_id: 2,
            in_nucleus: 1,
            dist_nucleus: 3.085635185241699,
            dist_periphery: 12.492626190185547,
        };
        for record in reader.records() {
            match record {
                Ok(r) => {
                    assert_eq!(r, expected_record);
                    break;
                }
                Err(_) => panic!(),
            }
        }
    }
}