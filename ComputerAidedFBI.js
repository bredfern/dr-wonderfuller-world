#!/usr/bin/nodejs


const x = require("./MyLib.js");

const util = require('util');
const exec = require('child_process').exec;


var xs = [];
for(let i = 0;i < 10;++i) xs.push(i);

for(let i = 0;i < 10;++i) {
	xs = x.shuffle(xs);
	console.log(xs);
}
	

let str = "";
process.argv.forEach((val, index) => {
	if(index > 1) str += " " + val + " ";
});

console.log("Eingabe ====> '" + str + "'");

let t = new x.tarot(str);
let r = new x.random(str);

function main() {
	//while(true) {
		let s = t.get();
		s = s.replace(/([0-9]+)[.]/gi,"$1ten");
		console.log("=====================================================");
		console.log(s);
		exec("trans -e google -b de:fr '" + s + "'",(err,stdout,stderr) => {
			if(stdout.trim() == "") {
				stdout = stdout.replace(/([0-9]+)[.]/gi,"$1ten");
				stdout = stdout.replace(/\\n/gi,"");
				console.log("Google => fail");
				exec("trans -e bing -b de:fr '" + s + "'",(err,stdout,stderr) => {
					stdout = stdout.replace(/([0-9]+)[.]/gi,"$1ten");
					stdout = stdout.replace(/\\n/gi,"")
					if(stdout.trim() == "") {
						console.log("Bing => fail");
						exec("/usr/bin/espeak -v mb-de7 -s 1 '" + s + "'",main);
					} else {
						console.log("Bing => success");
						exec("trans -e bing -b fr:de '" + stdout + "'",(err,stdout,stderr) => {
					if(stdout.trim() == "") {
						console.log("Bing back-translation => fail");
						exec("/usr/bin/espeak -v mb-de7 -s 1 '" + s + "'",main);
					} else {
							stdout = stdout.replace(/([0-9]+)[.]/gi,"$1ten");
							stdout = stdout.replace(/\\n/gi,"");
							console.log(stdout);
							exec("/usr/bin/espeak -v mb-de7 -s 1 '" + stdout + "'",main);
					}
						});
					}
				});
			} else { 
				console.log("Google => success");
				/*
				exec("trans -download-audio-as /tmp/xy -e google -b fr:de '" + stdout + "'",(err,stdout,stderr) => {
					console.log(stdout);
					exec("mpg123 /tmp/xy",main);
				});
				*/
				exec("trans -e google -b fr:de '" + stdout + "'",(err,stdout,stderr) => {
                                        if(stdout.trim() == "") {
                                                console.log("Google back-translation => fail");
                                                exec("/usr/bin/espeak -v mb-de7 -s 1 '" + s + "'",main);
                                        } else {
                                                        stdout = stdout.replace(/([0-9]+)[.]/gi,"$1ten");
                                                        stdout = stdout.replace(/\\n/gi,"");
                                                        console.log(stdout);
                                                        exec("/usr/bin/espeak -v mb-de7 -s 1 '" + stdout + "'",main);
                                        }
                               });

			}
		});
	//}
}


main();


