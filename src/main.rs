mod ising2d;

fn main() {
    let c = ising2d::solve_mc::<20>(10_000_000, 1.1);    
    println!("{:?}", c);
}
