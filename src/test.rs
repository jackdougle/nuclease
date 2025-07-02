pub fn run() {
    let mut lst: Vec<&str> = vec!["a", "b", "c"];
    lst.push("d");
    lst[2] = "King K. Rool";
    println!("{:?}", lst);

    for x in 0..lst.len() {
        let new: String = format!("{}: {}", x, lst[x]);
        println!("{}", new.as_str());
    }

    println!("{:?}", lst);
}
