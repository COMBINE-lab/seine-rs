use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::fs::*;
use std::io;
use std::io::prelude::*;
use std::io::{BufReader, Read};
use std::path::{Path, PathBuf};

use flate2::read::GzDecoder;

/*******************************************************************************/
/*                         Salmon Output Files                                 */
/*******************************************************************************/
#[derive(Debug)]
pub struct SalmonFiles {
    pub prefix: PathBuf,
    pub quant_file: PathBuf,
    pub bootstrap_file: PathBuf,
    pub ambig_file: PathBuf,
    pub mi_file: PathBuf,
    pub eq_file: PathBuf,
    pub names_tsv_file: PathBuf,
    pub cmd_file: PathBuf,
    pub collapsed_log_file: PathBuf,
    pub group_file: PathBuf,
    pub delta_file: PathBuf,
    pub cluster_file: PathBuf,
    pub gene_cluster_file: PathBuf,
}

// construct the files
impl SalmonFiles {
    pub fn new<P: AsRef<Path>>(dname: P) -> SalmonFiles {
        let dir = dname.as_ref();
        if !dir.exists() {
            panic!("The directory {} did not exist", dir.to_str().unwrap());
        }
        if !dir.is_dir() {
            panic!(
                "The path {} did not point to a valid directory",
                dir.to_str().unwrap()
            );
        }
        let aux_info = dir.join("aux_info");

        let mut eq_name = "eq_classes.txt";
        let mi_path = aux_info.join("meta_info.json");
        if mi_path.exists() {
            let file = File::open(mi_path);
            let reader = BufReader::new(file.unwrap());
            let jd: MetaInfo = serde_json::from_reader(reader).unwrap();

            eq_name = if jd.eq_class_properties.contains(&"gzipped".to_string()) {
                "eq_classes.txt.gz"
            } else {
                "eq_classes.txt"
            };
        }

        SalmonFiles {
            prefix: PathBuf::from(dir),
            ambig_file: aux_info.join("ambig_info.tsv"),
            mi_file: aux_info.join("meta_info.json"),
            quant_file: dir.join("quant.sf"),
            eq_file: aux_info.join(eq_name),
            bootstrap_file: aux_info.join("bootstrap").join("bootstraps.gz"),
            names_tsv_file: aux_info.join("bootstrap").join("names.tsv.gz"),
            cmd_file: dir.join("cmd_info.json"),
            cluster_file: dir.join("clusters.txt"),
            collapsed_log_file: dir.join("collapsed.log"),
            group_file: dir.join("groups.txt"),
            delta_file: dir.join("delta.log"),
            gene_cluster_file: dir.join("gene_cluster.log"),
        }
    }
}

/*******************************************************************************/
/*                         Equivalence Classes                                 */
/*******************************************************************************/

#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct MetaInfo {
    pub num_valid_targets: u32,
    pub serialized_eq_classes: bool,
    pub num_bootstraps: u32,
    pub num_eq_classes: u32,
    pub eq_class_properties: Vec<String>,
    pub samp_type: String,
}

#[derive(Debug)]
pub struct EqClass {
    pub labels: Vec<usize>,
    pub weights: Vec<f64>,
    // Todo: probably should be usize?
    pub count: u32,
}

#[derive(Debug)]
pub struct EqClassView<'a> {
    pub labels: &'a [usize],
    pub weights: &'a [f64],
    pub count: u32,
}

#[derive(Debug, Default)]
pub struct EqClassList {
    pub offsets: Vec<usize>,
    pub labels: Vec<usize>,
    pub weights: Vec<f64>,
    pub counts: Vec<u32>,
}

pub struct IterEqClassList<'a> {
    inner: &'a EqClassList,
    pos: usize,
}

/// Iterator over an EqList that yields EqClassView
impl<'a> Iterator for IterEqClassList<'a> {
    type Item = EqClassView<'a>;
    fn next(&mut self) -> Option<Self::Item> {
        self.pos += 1;
        self.inner.get(self.pos)
    }
}

impl EqClassList {
    pub fn iter(&self) -> IterEqClassList {
        IterEqClassList {
            inner: self,
            pos: 0,
        }
    }

    pub fn len(&self) -> usize {
        self.offsets.len() - 1
    }

    pub fn is_empty(&self) -> bool {
        self.offsets.len() == 1
    }

    pub fn push(&mut self, mut ec: EqClass) {
        let len = &ec.weights.len();
        self.offsets.push(self.offsets.last().unwrap() + len);
        self.labels.append(&mut ec.labels);
        self.weights.append(&mut ec.weights);
        self.counts.push(ec.count);
    }

    pub fn new() -> EqClassList {
        EqClassList {
            offsets: vec![0_usize],
            labels: Vec::<usize>::new(),
            weights: Vec::<f64>::new(),
            counts: Vec::<u32>::new(),
        }
    }

    pub fn get(&self, i: usize) -> Option<EqClassView> {
        if i + 1 >= self.offsets.len() {
            None
        } else {
            let p = self.offsets[i];
            let l = self.offsets[(i + 1)] - p;
            Some(EqClassView {
                labels: &self.labels[p..(p + l)],
                weights: &self.weights[p..(p + l)],
                count: self.counts[i],
            })
        }
    }
}

#[derive(Debug)]
pub struct EqClassCollection {
    pub targets: Vec<String>,
    pub ntarget: usize,
    pub neq: usize,
    pub classes: EqClassList,
}

impl EqClassCollection {
    // /// Add an equivalence class to the set of equivalence classes for this experiment
    // pub fn push(&mut self, ec: EqClass) {
    //     self.classes.push(ec);
    // }

    pub fn new() -> EqClassCollection {
        EqClassCollection {
            targets: Vec::<String>::new(),
            ntarget: 0,
            neq: 0,
            classes: EqClassList::new(),
        }
    }

    pub fn from_path<P: AsRef<Path>>(filename: &P) -> Result<EqClassCollection, io::Error> {
        let filename = filename.as_ref();
        let file = File::open(filename).expect("equivalence class file does not exist");
        let reader: Box<dyn Read> = if filename.ends_with("eq_classes.txt.gz") {
            Box::new(GzDecoder::new(file))
        } else {
            Box::new(file)
        };
        let mut buf_reader = BufReader::new(reader);
        let mut buf = String::new();

        let mut exp = EqClassCollection::new();

        buf_reader
            .read_line(&mut buf)
            .expect("Cannot read first line");
        buf.pop();
        let num_target: usize = buf.parse().unwrap();
        exp.ntarget = num_target;
        buf.clear();

        buf_reader
            .read_line(&mut buf)
            .expect("Cannot read second line");
        buf.pop();
        let num_eq: usize = buf.parse().unwrap();

        exp.neq = num_eq;

        let mut tnames = Vec::<String>::with_capacity(num_target);

        for _ in 0..num_target {
            buf.clear();
            buf_reader
                .read_line(&mut buf)
                .expect("could read target name");
            buf.pop();
            tnames.push(buf.to_string());
        }

        exp.targets = tnames;

        //let mut pb = pbr::ProgressBar::new(count);
        //pb.format("╢▌▌░╟");

        for _ in 0..num_eq {
            buf.clear();
            buf_reader
                .read_line(&mut buf)
                .expect("could read eq. class");
            buf.pop();
            let mut iter = buf.split_ascii_whitespace();
            let nt: usize = iter.next().unwrap().parse().unwrap();
            let mut tv = Vec::<usize>::with_capacity(nt);
            let mut wv = Vec::<f64>::with_capacity(nt);
            for _ in 0..nt {
                tv.push(iter.next().unwrap().parse().unwrap());
            }
            for _ in 0..nt {
                wv.push(iter.next().unwrap().parse().unwrap());
            }
            let c: u32 = iter.next().unwrap().parse().unwrap();

            let ec = EqClass {
                labels: tv,
                weights: wv,
                count: c,
            };
            exp.classes.push(ec);
            //pb.inc();
        }
        //pb.finish_print("done");
        Ok(exp)
    }

    pub fn get(&self, i: usize) -> Option<EqClassView> {
        self.classes.get(i)
    }
}

impl Default for EqClassCollection {
    fn default() -> Self {
        Self::new()
    }
}

/*******************************************************************************/
/*                         Quants                                              */
/*******************************************************************************/

#[derive(Debug, Deserialize, Serialize)]
pub struct QuantRecord {
    #[serde(rename = "Name")]
    pub name: String,

    #[serde(rename = "Length")]
    pub len: u32,

    #[serde(rename = "EffectiveLength")]
    pub efflen: f64,

    #[serde(rename = "TPM")]
    pub tpm: f64,

    #[serde(rename = "NumReads")]
    pub num_reads: f64,
}

#[derive(Debug, Deserialize)]
pub struct QuantEntry {
    pub len: u32,
    pub efflen: f64,
    pub tpm: f64,
    pub num_reads: f64,
}

/*******************************************************************************/
/*                         Extension Traits                                    */
/*******************************************************************************/

pub trait FromPathExt {
    fn from_path<P: AsRef<Path>>(p: P) -> Result<Self, csv::Error>
    where
        Self: Sized;
}

impl FromPathExt for HashMap<String, QuantEntry> {
    fn from_path<P: AsRef<Path>>(p: P) -> Result<HashMap<String, QuantEntry>, csv::Error> {
        // TODO error handling for csv::Error and io::Error
        let file = File::open(p);
        let mut rdr = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .from_reader(file.unwrap());

        let mut quant_map = HashMap::new();

        for quant_record in rdr.deserialize() {
            let quant_record: QuantRecord = quant_record?;
            let quant_entry = QuantEntry {
                len: quant_record.len,
                efflen: quant_record.efflen,
                tpm: quant_record.tpm,
                num_reads: quant_record.num_reads,
            };

            quant_map.insert(quant_record.name.to_string(), quant_entry);
        }

        Ok(quant_map)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ec_list_get() {
        let mut ecs = EqClassList::new();

        assert!(ecs.is_empty());

        let ec = EqClass {
            labels: vec![0, 1, 2],
            weights: vec![0.2, 0.3, 0.5],
            count: 15,
        };

        ecs.push(ec);

        assert_eq!(ecs.len(), 1);

        assert!(ecs.get(1).is_none());

        let ec = ecs.get(0).unwrap();
        assert_eq!(ec.labels, vec![0, 1, 2]);
        assert_eq!(ec.weights, vec![0.2, 0.3, 0.5]);
        assert_eq!(ec.count, 15);
    }
}
