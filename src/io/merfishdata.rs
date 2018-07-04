// Copyright 2016 Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

// TODO define trait for commonalities between tsv::Record and binary::Record
pub trait MerfishRecord {
    fn get_cell_id(&self) -> u32;
}

pub mod tsv {
    use std::io;
    use std::fs;
    use std::path::Path;

    use csv;
    use io::merfishdata::MerfishRecord;

    /// A 2D position in the microscope.
    #[derive(Serialize, Deserialize, Debug)]
    pub struct Position {
        #[serde(rename = "Position_X")]
        pub x: f32,
        #[serde(rename = "Position_Y")]
        pub y: f32,
    }


    /// A MERFISH raw data record.
    /// // "Cell_ID	Gene_Name	Hamming_Distance	Cell_Position_X	Cell_Position_Y	RNA_Position_X	RNA_Position_Y
    #[derive(Serialize, Deserialize, Debug)]
    pub struct Record {
        #[serde(rename = "Cell_ID")]
        pub cell_id: String,
        #[serde(rename = "Gene_Name")]
        pub feature: String,
        #[serde(rename = "Hamming_Distance")]
        pub hamming_dist: u8,
        #[serde(flatten, with = "prefix_cell")]
        pub cell_position: Position,
        #[serde(flatten, with = "prefix_rna")]
        pub rna_position: Position,
    }

    with_prefix!(prefix_cell "Cell_");
    with_prefix!(prefix_rna "RNA_");

    impl MerfishRecord for Record {
        fn get_cell_id(&self) -> u32 {
            self.cell_id.parse().expect("Failed parsing cell_id")
        }
    }


    /// A reader for MERFISH raw data.
    pub struct Reader<R: io::Read> {
        inner: csv::Reader<R>
    }


    impl Reader<fs::File> {
        /// Read from a given file path.
        pub fn from_file<P: AsRef<Path>>(path: P) -> io::Result<Self> {
            fs::File::open(path).map(Reader::from_reader)
        }
    }


    impl<R: io::Read> Reader<R> {
        pub fn from_reader(rdr: R) -> Self {
            Reader {
                inner: csv::ReaderBuilder::new().delimiter(b'\t').from_reader(rdr)
            }
        }

        pub fn records(&mut self) -> csv::DeserializeRecordsIter<R, Record> {
            self.inner.deserialize()
        }
    }
}

pub mod binary {
    use std::io;
    use std::io::Read;
    use std::fs;
    use bincode;
    use std::path::Path;
    use io::merfishdata::MerfishRecord;

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
        pub fn version(&self) -> u8 { self.version }
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
    #[derive(Serialize, Deserialize, Debug, PartialEq)]
    pub struct Record {
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

    impl MerfishRecord for Record {
        fn get_cell_id(&self) -> u32 {
            self.cell_id
        }
    }

    pub struct Reader<R: io::Read> {
        reader: io::BufReader<R>,
        header: Header,
        item: usize,
    }

    impl<R: io::Read> Reader<R> {
        /// Constructs a new `merfish::Reader<R>` where `R: io::Read`.
        ///
        /// This will read the header in order to check if the version matches.
        pub fn new(reader: R) -> Result<Self, &'static str> {
            let mut reader = io::BufReader::new(reader);
            let header: Result<Header, _> = bincode::deserialize_from(&mut reader);
            match header {
                Ok(ref header) if header.version != 1 => Err("Unsupported version {}"),
                Ok(ref header) if header.is_corrupt != 0 => Err("Header is corrupt, i.e. might not have been written properly (while in append mode)"),
                Ok(header) => {
                    // Read (discard) table header
                    let mut header_buf = vec![0u8; header.header_length as usize];
                    let r = reader.read_exact(&mut header_buf);
                    match r {
                        Err(_) => Err("Failed reading table header"),
                        Ok(()) => Ok(Reader {
                            reader,
                            header,
                            item: 0,
                        })
                    }
                }
                Err(_) => Err("Failed reading header")
            }
        }

        pub fn header(&self) -> &Header {
            &self.header
        }
    }

    impl Reader<fs::File> {
        /// Constructs a new `merfish::Reader<fs::File>` from given file.
        ///
        /// This delegates work to `merfish::Reader::new`.
        pub fn from_file<P: AsRef<Path>>(path: P) -> io::Result<Self> {
            fs::File::open(path).map(|r| { Reader::new(r).unwrap() })
        }
    }

    impl<R: io::Read> Iterator for Reader<R> {
        type Item = io::Result<Record>;

        /// Reads a single `merfish::Record` from the file, advancing the iterator.
        fn next(&mut self) -> Option<<Self as Iterator>::Item> {
            let (l, _) = self.size_hint();
            if l == 0 {
                None
            } else {
                match bincode::deserialize_from::<_, Record>(&mut self.reader) {
                    Ok(record) => {
                        self.item += 1;
                        Some(Ok(record))
                    }
                    // TODO: propagate bincode's Error
                    Err(_) => Some(Err(io::Error::new(io::ErrorKind::Other,
                                                      "Failed deserializing record from reader.")))
                }
            }
        }

        fn size_hint(&self) -> (usize, Option<usize>) {
            (self.header.num_entries as usize - self.item, Some(self.header.num_entries as usize - self.item))
        }
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use std::io;

    #[test]
    fn test_tsv_records() {
        let data = b"Cell_ID	Gene_Name	Hamming_Distance	Cell_Position_X	Cell_Position_Y	RNA_Position_X	RNA_Position_Y
0	SCUBE3	1	475.5	630.6	13146.86026973	25793.5656964
0	SCUBE3	1	475.5	630.6	13576.7356895	38396.4273422
";
        let mut reader = tsv::Reader::from_reader(io::Cursor::new(&data[..]));
        for r in reader.records() {
            match r {
                Ok(rec) => {
                    assert_eq!(rec.feature, "SCUBE3");
                    assert_eq!(rec.cell_position.x, 475.5);
                }
                Err(e) => panic!("{:?}", e)
            }
        }
    }

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

        let reader = binary::Reader::new(io::Cursor::new(&data[..])).unwrap();
        let expected_record = binary::Record {
            barcode: 4124,
            barcode_id: 1,
            fov_id: 0,
            total_magnitude: 22.515066,
            pixel_centroid: [174, 1513],
            weighted_pixel_centroid: [174.52765, 1512.7463],
            abs_position: [-1407.2377, -2946.629],
            area: 3,
            pixel_trace_mean: [0.0023688597, 0.13105613, 0.46080196, 0.2921718, 0.7618023, 0.007037689, 0.054068968, 0.0036673022, 0.000059783128, 0.0, 0.01174143, 0.0175779, 0.26575506, 0.107759856, 0.017091451, 0.07516251],
            pixel_trace_std: [0.0033500735, 0.04353166, 0.038341008, 0.02604062, 0.002876933, 0.0028763895, 0.049890008, 0.0028644174, 0.00008454611, 0.0, 0.0049137515, 0.024341747, 0.047965504, 0.04564901, 0.0126899, 0.037814822],
            is_exact: 1,
            error_bit: 0,
            error_dir: 0,
            av_distance: 0.46709287,
            cell_id: 2,
            in_nucleus: 1,
            dist_nucleus: 3.085635185241699,
            dist_periphery: 12.492626190185547,
        };
        for record in reader {
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
