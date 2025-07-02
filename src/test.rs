use rand::Rng;

pub fn run() {
    let chars: Vec<char> = vec!['A', 'T', 'C', 'G'];
    let mut rng = rand::rng;

    for x in 0..10 {
        let mut word = String::new();
        // Use rand::Rng's gen_range method and a proper rng instance
        while word.len() <= 120 {
            let n = rng().gen_range(0..4);
            word.push(chars[n]);
        }
        println!("{}", word);
    }
}
