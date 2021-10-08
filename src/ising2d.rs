use std::fmt;
use rand::{rngs::{StdRng}, RngCore, Rng, SeedableRng, distributions::{Distribution, Uniform}};

struct Ising2DContext<const L: usize> {
    adj : [[[(usize, usize); 4]; L]; L],
    range : Uniform<usize>,
    rand_rng : StdRng
}

impl<const L: usize> fmt::Debug for Ising2DContext<L> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {

        let mut result = String::from("Ising2DContext : \n{\n");
        for y in 0..L {
            for x in 0..L {
                result.push_str(&(format!("\t({}, {}) -> [ ", x, y)));
                for p in 0..4 {
                    let (a, b) = self.adj[x][y][p];
                    result.push_str(&(format!("({}, {}) ", a, b)));
                }
                result.push_str("]\n");
            }
        }
        result.push_str("}\n");

        return write!(f, "{}", result);
    }
}

fn initialize_context<const L: usize>() -> Ising2DContext<L> {
    let mut _adj = [[[(8888, 8888); 4]; L]; L];
    for (x, y) in itertools::iproduct!(0..L, 0..L) {
        _adj[x][y][0] = (((x + 1) % L), (y));
        _adj[x][y][1] = (((x + L - 1) % L), (y));
        _adj[x][y][2] = ((x), ((y + 1) % L));
        _adj[x][y][3] = ((x), ((y + L - 1) % L));
    }

    let mut seed = [128u8; 32];
    rand::thread_rng().fill_bytes(&mut seed);

    return Ising2DContext {
        adj : _adj,
        range : Uniform::from(0..L),
        rand_rng : StdRng::from_seed(seed)
    };
}

struct Ising2D<const L: usize> {
    spins : [[bool; L]; L],
    ham : f64
}

impl<const L: usize> fmt::Debug for Ising2D<L> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {

        let mut result = String::from("Ising2D : \n{\n");
        result.push_str("spins : \n\t{\n");
        for y in 0..L {
            result.push_str("\t\t[ ");
            for x in 0..L {
                result.push_str(if self.spins[x][y] { "+ " } else { "0 " });
            }
            result.push_str("]\n");
        }
        result.push_str("\t}\n");
        result.push_str(&(format!("ham : \t{}\n", self.ham)));

        return write!(f, "{}", result);
    }
}

fn interaction_energy<const L: usize>(context : &mut Ising2DContext<L>, spins : &[[bool; L]; L], (x, y) : (usize, usize), flip_spin : bool) -> f64 {
    let spin_energy = | s1 : bool, s2 : bool | {
        return if s1 == s2 { 1.0 } else { -1.0 };
    };
    
    let s1 = if flip_spin { !spins[x][y] } else { spins[x][y] };

    let mut h = 0.0;
    for n in 0..3 {
        let (nx, ny) = context.adj[x][y][n];
        let s2 = spins[nx][ny];
        h = h - spin_energy(s1, s2);
    }
    return h;
}

fn initialize_lattice<const L: usize>(context : &mut Ising2DContext<L>) -> Ising2D<L> {
    let mut spins = [[false; L]; L];

    for (x, y) in itertools::iproduct!(0..L, 0..L) {
        spins[x][y] = context.rand_rng.gen();
    }

    let ham = itertools::iproduct!(0..L, 0..L)
        .fold(0.0, | _h, _cell | {return _h + interaction_energy::<L>(context, &spins, _cell, false); });

    return Ising2D { spins, ham };
}

pub fn solve_mc<const L: usize>(num_iterations : i64, temp : f64) {
    let mut context = initialize_context::<L>();
    let mut lattice = initialize_lattice(&mut context);
    let mut flips = 0;

    let evolve = | context: &mut Ising2DContext<L>, lattice: &mut Ising2D<L>, f: &mut i64, cell: (usize, usize) | {
        let before = interaction_energy::<L>(context, &lattice.spins, cell, false);
        let after  = interaction_energy::<L>(context, &lattice.spins, cell, true);
        let d_e = after - before;
        let accept = 
            if d_e <= 0.0 { 
                true
            } else {
                *f += 1;
                let probability = (-d_e/temp).exp();
                let random = context.rand_rng.gen::<f64>();
                probability >= random
            };

        if accept {
            lattice.ham = lattice.ham + d_e;
            let (x, y) = cell;
            lattice.spins[x][y] = !lattice.spins[x][y];
        }
    };
    
    let start = std::time::Instant::now();
    println!("Before: {:?}", lattice);
    for _ in 0..num_iterations {
        let x = context.range.sample(&mut context.rand_rng);
        let y = context.range.sample(&mut context.rand_rng);
        evolve(&mut context, &mut lattice, &mut flips, (x, y));
    }
    let finish = std::time::Instant::now();
    println!("After: {:?}", lattice);
    
    let duration = finish.duration_since(start);
    println!("Rust - Solving Ising2D {side} x {side} ({num_iterations} iterations) at {temp}K took {duration}ms. {flips} spins were flipped.", side=L, num_iterations=num_iterations, temp=temp, duration=duration.as_millis(), flips=flips);
}