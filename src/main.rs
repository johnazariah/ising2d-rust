mod ising2d;

fn main() {
    let c = ising2d::solve_mc::<70>(5_000_000, 0.6);    
    println!("{:?}", c);
}
