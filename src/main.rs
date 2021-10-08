mod ising2d;

fn main() {
    let c = ising2d::solve_mc::<70>(250_000_000, 1.2);    
    println!("{:?}", c);
}
