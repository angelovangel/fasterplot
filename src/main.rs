use std::{collections::BTreeMap, process::exit};
use clap::{Parser, Args};
use clap::builder::styling::*;
use kseq::parse_path;
extern crate histo;
use histo::Histogram;

// here at least one of Flags is required, multiple are not allowed, and infile is positional required
#[derive(Parser)]
#[command(
    version, about, long_about = None, 
    // clap v4 is not colored anymore, have to be explicit
    styles = Styles::styled()
        .header(AnsiColor::Yellow.on_default())
        .usage(AnsiColor::Yellow.on_default())
        .literal(AnsiColor::Green.on_default())
        .placeholder(AnsiColor::Green.on_default())
)
]
struct Cli {
    #[arg(required = true)]
    infile: Option<String>,
    
    #[command(flatten)]
    flags: Flags,
    
    #[arg(short ='x', long, required = false, help = "output histogram instead of raw data")]
    hist: bool,
    
    #[arg(short ='s', 
        long = "skip", 
        required = false, 
        default_value_t = 1,
        value_parser = clap::value_parser!(u16).range(0..1001),
        help = "skip every nth read [0..1000] to speedup processing for large files")
    ]
    skip: u16,

    #[arg(short, long, 
        required = false,
        default_value_t = 50000, 
        value_parser = clap::value_parser!(u32).range(100..500001),
        help = "max length to use in len output [1000..500000]")
    ]
    maxlen: u32,
    
    #[arg(short, long, 
        required = false, 
        default_value_t = 500,
        value_parser = clap::value_parser!(u32).range(10..1001),
        help = "bin step in len output [10..1000]")
    ]
    window: u32,
}

#[derive(Args)]
#[group(required = true, multiple = false)]
struct Flags {
    #[arg(short, long, help = "output base yield over qscores")]
    qscore: bool,
    #[arg(short, long, help = "output base yield over length bins")]
    len: bool, 
    #[arg(short, long, help = "output Nx values (from N10 to N100)")]
    nx: bool,
}


// struct to hold read counts and bases per len bin
// struct Lenbin {
//     readcount: i64,
//     bases: i64,
// } 
fn main() {
    let cli = Cli::parse();
    let fastqfile = cli.infile.unwrap();
    let mut records = parse_path(fastqfile).unwrap();
    let mut bases: i64 = 0;

    // base yield over qscores
    if cli.flags.qscore {
        let mut hash: BTreeMap<u8, i64> = BTreeMap::new();
        let qvals: Vec<u8> = (33..=126).step_by(1).collect();
        //let vals: Vec<u8> = (33..=83).step_by(1).collect(); //0 to 50 only
        for v in qvals {
            hash.insert(v, 0);
        }

        // hist all qvals for all bases
        let mut histogram = Histogram::with_buckets(10);
        let mut step: u16 = 0;
        while let Some(record) = records.iter_record().unwrap() {
            step += 1;
            if step != cli.skip {
                continue;
            }
            step = 0;
            bases += record.len() as i64;
            let quals = record.qual().as_bytes().to_owned();
            for q in quals {
                if cli.hist{
                    histogram.add(q as u64 - 33)
                }
                *hash.entry(q).or_insert(0) += 1;
                //hash.insert(q,  hash[&q] + 1);
            }
        }
        if cli.hist{
            println!("# Histogram of qvalues for all bases");
            println!("{}", histogram);
            exit(0)
        }

        // print qyields
        println!("qvalue\tbases_at_q\tbases_above_q\tpercent_at_q\tpercent_above_q");
        let mut cumbases: i64 = bases;
        for (key, value) in hash.range(33..=83) { //0 to 50 only
            cumbases -= value;
            println!(
                "{}\t{}\t{}\t{}\t{}", 
                key - 33, value, cumbases, 
                format!("{:.4}", *value as f64 / bases as f64 * 100.0),
                format!("{:.4}", cumbases as f64 / bases as f64 * 100.0)
            );
        }
    }
    // base yield over lengths
    if cli.flags.len {
        let mut lens:Vec<i32> = Vec::new();
        let mut step: u16 = 0;

        while let Some(record) = records.iter_record().unwrap() {
            step += 1;
            if step != cli.skip {
                continue;
            }
            step = 0;
            bases += record.len() as i64;
            lens.push(record.len() as i32);
            //lens.sort_unstable();
        } 
        // do not allow cli.maxlen to be smaller than read length
        let maxlen = *lens.iter().max().unwrap() as u32;
        if cli.maxlen < maxlen + cli.window {
            println!(
                "Max read length is {}, with --window {} please choose --maxlen bigger than {}", 
                maxlen, cli.window, maxlen + cli.window
            );
            exit(1)
        }

        if cli.hist {
            let mut histogram = Histogram::with_buckets(10);
            for l in lens {
                histogram.add(l as u64)
            }
            println!("# Histogram of read lengths");
            println!("{}", histogram);
            exit(0);
        }
        let mut lenhash: BTreeMap<i32, i64> = BTreeMap::new();
        let lenvals:Vec<u32> = (1..cli.maxlen).step_by(cli.window as usize).collect();
        let mut cumbases = bases;
        for k in lenvals {
            lenhash.insert(k as i32, 0);
        }

        for l in lens{
            // get the next back key for this read length
            if l <= 1 { // crashes if l == 1, so 
                //lenhash.insert(1, 1);
                continue;
            }
            let mykey = lenhash.range(l..).next().unwrap().0;
            *lenhash.entry(*mykey).or_insert(0) += l as i64;
        }
        
        // print output
        let maxbin_bases = lenhash.iter().max_by_key( |p| p.1 ).unwrap().1; 
        println!("# maxbin:\t{}", lenhash.iter().max_by_key( |p| p.1 ).unwrap().0);
        println!("# total_bases:\t{}", bases);
        println!("# maxbin_bases:\t{}", maxbin_bases);
        println!("lenbin\tbases\tbases_above_len\tpercent_at_lenbin\tpercent_above_lenbin");
        for (key, value) in lenhash {
            cumbases -= value;
            println!(
                "{}\t{}\t{}\t{}\t{}", 
                key, value, cumbases, 
                format!("{:.4}", value as f64 / bases as f64 * 100.0),
                format!("{:.4}", cumbases as f64 / bases as f64 * 100.0)
            );
        }
    }
    
    // Nx values
    if cli.flags.nx {
        let nbins:Vec<i16> = (0..10).step_by(1).collect();
        let mut lens:Vec<i64> = Vec::new();

        while let Some(record) = records.iter_record().unwrap() {
            bases += record.len() as i64;
            lens.push(record.len() as i64);
            
        }

        lens.sort_unstable();
        let cumsum = lens
            .iter()
            .scan(0, |sum, i|{
                *sum += i;
                Some(*sum)
            }).collect::<Vec<_>>();
        
        // output
        println!("Nx\tread_len");
        for i in nbins {
            let nx = i as f32 / 10.0;
            let nxsum = lens.iter().sum::<i64>() as f32 * nx;
            let nxindex = cumsum.iter().position(|&x| x > nxsum as i64).unwrap();
            println!("{}\t{}", format!("{:.0}", 100.0 - nx * 100.0), lens[nxindex]);

        }
    }

}