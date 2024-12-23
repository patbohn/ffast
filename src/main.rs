use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::PathBuf;
use flate2::read::MultiGzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;

use structopt::StructOpt;

///get phred (intergers) to error probability table
pub mod phred_int_to_prob;
use phred_int_to_prob::PHRED_TO_ERROR_PROB;

#[derive(Debug, StructOpt)]
#[structopt(name = "calc_per_read_rust", about = "Calculates the average quality score and error probability per base per read in a gzipped fastq file.")]
struct Config {
    #[structopt(parse(from_os_str), help = "Input fastq file")]
    input_file: PathBuf,

    #[structopt( help = "Output csv file with read ID, length, mean quality score, and error probability per base", default_value ="output")]
    output_prefix: String,
}

fn main() -> std::io::Result<()> {
    let config = Config::from_args();
    
    let per_read_filename = format!("{}_read_stats.csv.gz", config.output_prefix);
    let hist_filename = format!("{}_phred_hist.csv", config.output_prefix);
    
    let input_file = File::open(config.input_file)?;
    let reader = BufReader::new(MultiGzDecoder::new(input_file));
    
    let per_read_output_file = File::create(per_read_filename)?;
    
    let mut writer = GzEncoder::new(per_read_output_file, Compression::fast());

    writeln!(&mut writer, "read_id,read_length,mean_phred,mean_error_rate")?;
    
    let mut lines_iter = reader.lines();
    
    //create array to count phred qscores
    let mut qscore_hist: [u64; 100] = [0;100];
    
    
    while let Some(Ok(header)) = lines_iter.next() {
        let _sequence = lines_iter.next().unwrap().unwrap();
        let _ = lines_iter.next(); // ignore separator line
        let quality = lines_iter.next().unwrap().unwrap();

        for c in quality.chars().filter(|&c| c >= '\x21' && c <= '\x7f') {
            
            qscore_hist[c as usize - 33] += 1;
        }
        
        let id = get_read_id(&header);
        let (mean_error_prob, read_length) = calc_mean_median_error(&quality);
        let mean_quality = error_prob_to_phred(mean_error_prob);
        
        writeln!(&mut writer, "{},{}, {:.2},{:1.2e}", id, read_length, mean_quality, mean_error_prob)?;
    }
    
    //write phred qscore histogram
    let hist_output_file = File::create(hist_filename)?;
    let mut writer = BufWriter::new(hist_output_file);
    writeln!(&mut writer, "phred_score,count")?;

    for (qscore, count) in qscore_hist.iter().enumerate() {
        writeln!(&mut writer, "{},{}", qscore, count)?;
        }

    Ok(())
}


fn get_read_id(read_id_line: &str) -> &str {
    match read_id_line.split_once(" ") {
        Some((read_id, _additional_info)) => {
            return read_id.chars().next().map(|c| &read_id[c.len_utf8()..]).unwrap_or("NaN");
            }
        None => {return read_id_line}
        }
}

fn calc_mean_median_error(quality_str: &str) -> (f64, i64) {
    let mut total_prob = 0.0;
    let mut count = 0;
    
    for quality_char in quality_str.chars() {
        let phred = (quality_char as i8) - 33;
        let prob = PHRED_TO_ERROR_PROB[phred as usize];
        total_prob += prob;

        
        count += 1;
    }

    
    return (total_prob / count as f64, count as i64);
}


//replaced by lookup table
/*
fn phred_to_error_prob(phred: f64) -> f64 {
    10.0_f64.powf(-phred / 10.0)
}
*/

fn error_prob_to_phred(prob: f64) -> f64 {
    return -10.0_f64 * prob.log10()
}

