

let kernel = null;
if(typeof gpu !== "undefined") kernel = gpu.createKernel(instrumentF).setOutput([50000]);









function fractal(cfg) {
	let mathx = (cfg[3] - cfg[2]) * this.thread.x / cfg[0] + cfg[2];
	let mathy = cfg[5] - (cfg[5] - cfg[4]) * this.thread.y / cfg[1];
	let xr = 0;
	let xi = 0;
	for(let k = 0;k < 100;++k) {
		let zr = 1;
		let zi = 0;
		
		let counterPolyR = mathx;
		let counterPolyI = mathy;
		for(let i = 0;i < 5;++i) {
			let zr2 = zr * xr - zi * xi;
			zi = zr * xi + zi * xr;
			zr = zr2;
			
			counterPolyR += cfg[6 + 2 * i] * zr - cfg[6 + 2 * i + 1] * zi;
			counterPolyI += cfg[6 + 2 * i] * zi + cfg[6 + 2 * i + 1] * zr;
		}

		zr = 1;
		zi = 0;
		
		let denominatorPolyR = mathx;
		let denominatorPolyI = mathy;
		
		for(let i = 0;i < 5;++i) {
			let zr2 = zr * xr - zi * xi;
			zi = zr * xi + zi * xr;
			zr = zr2;
			
			denominatorPolyR += cfg[16 + 2 * i] * zr - cfg[16 + 2 * i + 1] * zi;
			denominatorPolyI += cfg[16 + 2 * i] * zi + cfg[16 + 2 * i + 1] * zr;
		}
		
		let z = denominatorPolyR * denominatorPolyR + denominatorPolyI * denominatorPolyI;
		xr = (counterPolyR * denominatorPolyR + counterPolyI * denominatorPolyI) / z;
		xi = (-counterPolyR * denominatorPolyI + counterPolyI * denominatorPolyR) / z;

		if(xr * xr + xi * xi > 10.0) {
			this.color(
				(Math.sin(0.05 * k * 2.0 * Math.PI) + 1.0) / 2.0,
				(Math.sin(0.06 * k * 2.0 * Math.PI) + 1.0) / 2.0,
				(Math.sin(0.07 * k * 2.0 * Math.PI) + 1.0) / 2.0,
				1
			);
			return k;
		}
	}
	this.color(0,0,0,0);
}



function shuffle(xs) {
xs = xs.slice(0);

	var ys = [];
	
	while(xs.length > 0) {
		var i = Math.round(Math.random() * (xs.length - 1));
		ys.push(xs[i]);
		xs.splice(i,1);
	}
	
	return ys;
}


class Chess {
	team() {
		return this._team;
	}

	write() {
		document.write("<legend>" + this.name() +
			 "  (<input type='checkbox' id='" + this.cvID + "3'> Coop)  "
			 + "</legend>");
		document.write("<canvas width=300 height=300 id='" + this.cvID + "'></canvas>");
		document.write("<br><table><tr><td><select id='" + this.cvID + "1'><option value=0>0</option><option value=1>1</option><option value=2>2</option><option value=3>3</option><option value=4>4</option></select>");
		document.write("<button id='" + this.cvID + "2'>Vorschlag</button>");
		document.write("</td><td id='" + this.cvID + "0'></td></tr></table>");
		
		this.render();
		
		document.getElementById(this.cvID + "3").addEventListener("click",((evt) => {
			if(!this.thinking)
				this._team = document.getElementById(this.cvID + "3").checked;
			document.getElementById(this.cvID + "3").checked = this._team;
		}).bind(this));

		document.getElementById(this.cvID + "2").addEventListener("click",((evt) => {
			if(this.thinking) return;
			this.propose(parseInt(document.getElementById(this.cvID + '1').value));
		}).bind(this));
		
		document.getElementById(this.cvID).addEventListener("mousemove",((evt) => {
			if(this.thinking) return;
			/*if(onMove) return;
			if(isMouseDown) return;*/
			var cv = document.getElementById(this.cvID);
			var xy = this.getOffset(cv);
			
			this.mousex = Math.floor((evt.x - xy.left) * 8 / cv.width);
			this.mousey = Math.floor((evt.y - xy.top) * 8 / cv.height);

			if(this.mousex > 7) this.mousex = 7;
					
			if(!this.board[this.mousex + 8 * this.mousey] || (this.board[this.mousex + 8 * this.mousey].color != this.board.color && !this.team())) {
				this.mousex = -1;
				this.mousey = -1;
			}
			
			this.render();
		}).bind(this));

		document.getElementById(this.cvID).addEventListener("click",((evt) => {
			if(this.thinking) return;
			/*isMouseDown = false;
			if(onMove) return;*/
			if(this.measureWhite)
				if(document.getElementById(this.cvID + "0"))
					document.getElementById(this.cvID + "0").innerHTML = "turn(" + this.board.color + ")  <br> points(" + this.measureWhite(this.board) + ")";

			if(this.mousex < 0 && this.selectedx < 0) return;		
			
			var cv = document.getElementById(this.cvID);
			var xy = this.getOffset(cv);
			
			let sx = Math.floor((evt.x - xy.left) * 8 / cv.width);
			let sy = Math.floor((evt.y - xy.top) * 8 / cv.height);

			if(sx == this.selectedx && sy == this.selectedy) {
				this.selectedMoves = [];
				this.selectedx = -1;
				this.selectedy = -1;
				this.render();
				return;
			}
			
			if(this.selectedMoves.length == 1)
			for(var i = 0;i < this.selectedMoves[0][2].length;++i) {
				var m = this.selectedMoves[0][2][i];
				
				if(sx == m[0] && sy == m[1]) {
					this.board = this.doMove(this.board,this.selectedx,this.selectedy,m);
					this.selectedx = -1;
					this.selectedy = -1;
					this.mousex = -1;
					this.mousey = -1;
					this.selectedMoves = [];
					this.lastMoveInitialX = this._lastMoveInitialX;
					this.lastMoveInitialY = this._lastMoveInitialY;
					this.lastMoveX = sx;
					this.lastMoveY = sy;
					this.pp = null;
					this.pp2 = null;
					this.render();
					if(this.measureWhite)
						if(document.getElementById(this.cvID + "0"))
							document.getElementById(this.cvID + "0").innerHTML = "turn(" + this.board.color + ")  <br> points(" + this.measureWhite(this.board) + ")";
					return;
				}
			}
			
			this.selectedx = sx;
			this.selectedy = sy;
		
			if(!this.board[this.selectedx + 8 * this.selectedy]) {
				this.selectedx = -1;
				this.selectedy = -1;
				this.selectedMoves = [];
			} else if(this.team() || this.board[this.selectedx + 8 * this.selectedy].color == this.board.color){
				this._lastMoveInitialX = this.selectedx;
				this._lastMoveInitialY = this.selectedy;
				this.selectedMoves = [[sx,sy,this.endFilter(this.board,1,this.selectedx,this.selectedy,this.moveFilter(this.board,this.selectedx,this.selectedy,this.getMoves(this.board,this.selectedx,this.selectedy)))]];
			} else {
				this.selectedx = -1;
				this.selectedy = -1;
				this.selectedMoves = [];
			}
			
			this.render();
		}).bind(this));

		if(this.measureWhite)
			if(document.getElementById(this.cvID + "0"))
				document.getElementById(this.cvID + "0").innerHTML = "turn(" + this.board.color + ")  <br> points(" + this.measureWhite(this.board) + ")";

	}
	
	constructor() {
		this.cvID = "" + Math.random();
		this.ps = {};
		this.ps["queen"] = 1000;
		this.ps["rook"] = 100;
		this.ps["knight"] = 10;
		this.ps["bishop"] = 10;
		this.ps["pawn"] = 1;
		this.ps["king"] = 10000;
		
		this.init(this.cvID);
		
		this.selectedMoves = [];
		this.mousex = -1;
		this.mousey = -1;
		this.onMove = true;
		this.lastMoveX = -1;
		this.lastMoveY = -1;
		this.selectedx = -1;
		this.selectedy = -1;
		this.lastMoveInitialX = -1;
		this.lastMoveInitialY = -1;
		
		this.board.color = "white";
		this.board.counter = 0;
	}
	
	getOffset( el ) {
		var _x = 0;
		var _y = 0;
		while( el && !isNaN( el.offsetLeft ) && !isNaN( el.offsetTop ) ) {
		_x += el.offsetLeft - el.scrollLeft;
		_y += el.offsetTop - el.scrollTop;
		el = el.offsetParent;
		}
		return { top: _y, left: _x };
	}

	loadPiece(name,imgName,color) {
		var y = {};
		
		if(typeof document !== "undefined") y = document.createElement("img");
		y.src = "img/" + imgName + ".png";
		y.shape = name;
		y.color = color;
		y.onload = this.render.bind(this);
		
		return y;
	}

	render() {
		var cv = document.getElementById(this.cvID);
		if(cv == null) return;
		var ctx = cv.getContext("2d");
		
		ctx.fillStyle = "black";
		ctx.fillRect(0,0,cv.width,cv.height);
		
		var c = 0;
		for(var i = 0;i < 8;++i) {
			for(var j = 0;j < 8;++j) {
				ctx.fillStyle = c ? "lightgray" : "white";
				
				/*for(var k = 0;k < this.selectedMoves.length;++k)
					if(i == this.selectedMoves[k][0] && j == this.selectedMoves[k][1]) {
						ctx.fillStyle = "violet";
					}*/
					
				ctx.fillRect(i * cv.width / 8,j * cv.height / 8,cv.width / 8,cv.height / 8);

					
				c = 1 - c;
			}
			c = 1 - c;
		}
		
		for(var i = 0;i < 8;++i) {
			for(var j = 0;j < 8;++j) {
				if(this.board[i + 8 * j] != null) 
					ctx.drawImage(this.board[i + 8 * j],i * cv.width / 8,j * cv.height / 8,cv.width / 8,cv.height / 8);
					if(this.pp2)
					if(this.pp2[i + "," + j]) {
						ctx.strokeStyle = "magenta"
						ctx.lineWidth = 4;
						
						ctx.beginPath();
						ctx.rect(i * cv.width / 8,j * cv.height / 8,cv.width / 8,cv.height / 8);
						ctx.stroke();
					}
			}
		}
		
		/*if(this.pp)
		Object.keys(this.pp[1]).forEach(p => {
			let ps = p.split(",");
			
			ctx.strokeStyle = "cyan"
			ctx.lineWidth = 5;
			ctx.beginPath();
			ctx.rect(parseInt(ps[0]) * cv.width / 8,parseInt(ps[1]) * cv.height / 8,cv.width / 8,cv.height / 8);
			ctx.stroke();
			
			ctx.strokeStyle = "magenta"
			ctx.lineWidth = 5;
			ctx.beginPath();
			ctx.rect(parseInt(ps[2]) * cv.width / 8,parseInt(ps[3]) * cv.height / 8,cv.width / 8,cv.height / 8);
			ctx.stroke();
		});*/
		
		if(this.lastMoveX >= 0) {
			ctx.strokeStyle = "yellow"
			ctx.lineWidth = 2;
			
			ctx.beginPath();
			ctx.rect(this.lastMoveX * cv.width / 8,this.lastMoveY * cv.height / 8,cv.width / 8,cv.height / 8);
			ctx.stroke();
		}

		if(this.lastMoveInitialX >= 0) {
			ctx.strokeStyle = "brown"
			ctx.lineWidth = 2;
			
			ctx.beginPath();
			ctx.rect(this.lastMoveInitialX * cv.width / 8,this.lastMoveInitialY * cv.height / 8,cv.width / 8,cv.height / 8);
			ctx.stroke();
		}

		for(var i = 0;i < this.selectedMoves.length;++i) {
			for(var k = 0;k < this.selectedMoves[i][2].length;++k) {
				var x = this.selectedMoves[i][2][k];
				
				ctx.strokeStyle = "blue"
				ctx.lineWidth = 6;
				
				ctx.beginPath();
				ctx.rect(x[0] * cv.width / 8,x[1] * cv.height / 8,cv.width / 8,cv.height / 8);
				ctx.stroke();
			}
		}
		
		for(var i = 0;i < this.selectedMoves.length;++i) {
			for(var k = 0;k < this.selectedMoves[i][2].length;++k) {
				var x = this.selectedMoves[i][2][k];
				if(this.pp)
				if(this.pp[1][
					this.selectedMoves[i][0] + "," + this.selectedMoves[i][1]
					+ "," + x[0] + "," + x[1]]
				) {
					ctx.strokeStyle = "cyan"
					ctx.lineWidth = 2;
					
					ctx.beginPath();
					ctx.rect(x[0] * cv.width / 8,x[1] * cv.height / 8,cv.width / 8,cv.height / 8);
					ctx.stroke();
				}
			}
		}
		
		if(this.mousex >= 0) {
			ctx.strokeStyle = "green"
			ctx.lineWidth = 5;
			
			ctx.beginPath();
			ctx.rect(this.mousex * cv.width / 8,this.mousey * cv.height / 8,cv.width / 8,cv.height / 8);
			ctx.stroke();
		}
		
		if(this.selectedx >= 0) {
			ctx.strokeStyle = "red"
			ctx.lineWidth = 5;
			
			ctx.beginPath();
			ctx.rect(this.selectedx * cv.width / 8,this.selectedy * cv.height / 8,cv.width / 8,cv.height / 8);
			ctx.stroke();
		}
		
		if(this.progress != null) {
			ctx.strokeStyle = "blue"
			ctx.lineWidth = 5;
			
			ctx.beginPath();
			ctx.rect(1 * cv.width / 8,3.5 * cv.height / 8,6 * cv.width / 8,1 * cv.height / 8);
			ctx.stroke();

			ctx.fillStyle = "blue"
			ctx.fillRect(1 * cv.width / 8,3.65 * cv.height / 8,this.progress * 6 * cv.width / 8 / 100,0.7 * cv.height / 8);
		}
		
		/*
		if(generalInit == gardenEdenInit) {
			var c = this.board[0].color;
			var cnt = 0;
			
			for(var i = 0;i < 8;++i) {
				if(this.board[i].color == c) cnt += 1;
				
				c = c == "white" ? "black" : "white";
			}
			for(var i = 0;i < 6;++i) {
				if(this.board[15 + 8 * i].color == c) cnt += 1;
				
				c = c == "white" ? "black" : "white";
			}
			for(var i = 0;i < 8;++i) {
				if(this.board[63 - i].color == c) cnt += 1;
				
				c = c == "white" ? "black" : "white";
			}
			for(var i = 0;i < 6;++i) {
				if(this.board[48 - 8 * i].color == c) cnt += 1;
				
				c = c == "white" ? "black" : "white";
			}
			
			ctx.fillStyle = "blue";
			ctx.font = "48px serif";
			ctx.fillText(Math.round(100 * cnt / 28) + "% complete",cv.width / 4,cv.height / 3);
		}*/
	}

	
	getMoves(board,i,j) {
		var piece = board[i + 8 * j];
		
		var xs = [];
		
		if(piece != null)
		switch(piece.shape) {
			case "pawn":
				if(piece.color == "black") {
					if(!board[i + 8 * (j + 1)]) xs.push([i,j + 2]);
					xs.push([i + 1,j + 1]);		
					xs.push([i - 1,j + 1]);		
					xs.push([i + 1,j - 1]);		
					xs.push([i - 1,j - 1]);		
					xs.push([i + 0,j + 1]);		
					xs.push([i + 0,j - 1]);		
					xs.push([i + 1,j + 0]);		
					xs.push([i - 1,j - 0]);		
				} else {
					if(!board[i + 8 * (j - 1)]) xs.push([i,j - 2]);
					xs.push([i + 1,j + 1]);		
					xs.push([i - 1,j + 1]);		
					xs.push([i + 1,j - 1]);		
					xs.push([i - 1,j - 1]);		
					xs.push([i + 0,j + 1]);		
					xs.push([i + 0,j - 1]);		
					xs.push([i + 1,j + 0]);		
					xs.push([i - 1,j - 0]);		
				}
				break;
			case "rook":
				for(var k = 0;k < 8;++k) {
					xs.push([i,j + k + 1]);
					xs.push([i,j - k - 1]);
					xs.push([i + k + 1,j]);
					xs.push([i - k - 1,j]);
				}
				break;
			case "bishop":
				for(var k = 0;k < 8;++k) {
					xs.push([i + k + 1,j + k + 1]);
					xs.push([i - k - 1,j - k - 1]);
					xs.push([i - k - 1,j + k + 1]);
					xs.push([i + k + 1,j - k - 1]);
				}
				break;
			case "queen":
				for(var k = 0;k < 8;++k) {
					xs.push([i,j + k + 1]);
					xs.push([i,j - k - 1]);
					xs.push([i + k + 1,j]);
					xs.push([i - k - 1,j]);
					xs.push([i + k + 1,j + k + 1]);
					xs.push([i - k - 1,j - k - 1]);
					xs.push([i - k - 1,j + k + 1]);
					xs.push([i + k + 1,j - k - 1]);
				}
				break;	
			case "knight":
				xs.push([i + 1,j + 2]);		
				xs.push([i - 1,j + 2]);		
				xs.push([i + 1,j - 2]);		
				xs.push([i - 1,j - 2]);		
				xs.push([i + 2,j + 1]);		
				xs.push([i - 2,j + 1]);		
				xs.push([i + 2,j - 1]);		
				xs.push([i - 2,j - 1]);	
				break;
					
			case "king":
				xs.push([i + 1,j + 1]);		
				xs.push([i - 1,j + 1]);		
				xs.push([i + 1,j - 1]);		
				xs.push([i - 1,j - 1]);		
				xs.push([i + 0,j + 1]);		
				xs.push([i + 0,j - 1]);		
				xs.push([i + 1,j + 0]);		
				xs.push([i - 1,j - 0]);		
				
				var brd = board;
				
				if(
					i == 4 && j == 0 && brd[4] && brd[4].color == "black" && brd[4].shape == "king" 
					&& !brd[5] && !brd[6] && brd[7] && brd[7].color == "black" && brd[7].shape == "rook"
				) {
					xs.push([6,0]);
				}
				
				if(
					i == 4 && j == 0 && brd[4] && brd[4].color == "black" && brd[4].shape == "king" 
					&& !brd[3] && !brd[2] && !brd[1] && brd[0] && brd[0].color == "black" && brd[0].shape == "rook"
				) {
					xs.push([2,0]);
				}
				
				if(
					i == 4 && j == 7 && brd[60] && brd[60].color == "white" && brd[60].shape == "king" 
					&& !brd[61] && !brd[62] && brd[63] && brd[63].color == "white" && brd[63].shape == "rook"
				) {
					xs.push([6,7]);
				}
				
				if(
					i == 4 && j == 7 && brd[60] && brd[60].color == "white" && brd[60].shape == "king" 
					&& !brd[59] && !brd[58] && !brd[57] && brd[56] && brd[56].color == "white" && brd[56].shape == "rook"
				) {
					xs.push([2,7]);
				}
		}
		
		return xs;
	}


	isUnderThreat(board,n,i,j) {
		for(var a = 0;a < 8;++a) {
			for(var b = 0;b < 8;++b) {
				if(board[a + 8 * b] && (board[a + 8 * b].color != board[i + 8 * j].color)) {
					var ms = this.moveFilter(board,a,b,this.getMoves(board,a,b));
					
					for(var k = 0;k < ms.length;++k) {
						if((ms[k][0] == i) && (ms[k][1] == j))
							return true;
					}
				}
			}
		}
		
		return false;
	}


	shuffle(xs) {
		var ys = [];
		
		while(xs.length > 0) {
			var i = Math.round(Math.random() * (xs.length - 1));
			ys.push(xs[i]);
			xs.splice(i,1);
		}
		
		return ys;
	}
	
	next(board) {
		if(board.color == "white")
			return "black";
		else
			return "white";
	}
	
	doMove(board,i,j,m) {
		let b = this._doMove(board,i,j,m);
		b.color = board.color;
		b.color = this.next(b);
		b.counter = board.counter + 1;
		
		return b;
	}
	
	cloneBoard(b) {
		let b2 = [];
		
		b.forEach(p => {
			if(p)
				b2.push({
					color : p.color,
					shape : p.shape
				})
			else
				b2.push(null);
		});
		
		b2.color = b.color;
		b2.counter = b.counter;
		
		return b2;
	}
	
	async _propose(board,i0,j0,mv0,mv1,n,zs) {
		if(n <= 2) {
			let xy = window.setInterval((() => {
				for(let k = 0;k < this.workers.length;++k) {
					let w = this.workers[k];
					if(w.done) {
						window.clearInterval(xy);
						w.done = false;
						w.onmessage = ((m) => {
							let xs = {};
							xs[[i0,j0,mv0,mv1] + ""] = 1;
							zs.push([m.data[1],xs]);
							w.done = true;
						}).bind(this);
						w.postMessage({
							type : this.constructor.name,
							thisBoard : this.cloneBoard(this.board),
							n : n,
							_team : this.team(),
							board : this.cloneBoard(board),
						});
						break;
					}
				}
			}).bind(this),100);
			
			return;
		}

		let num = 0;
		let zs2 = [];

		for(let i = 0;i < 8;++i)
		for(let j = 0;j < 8;++j)
		if(board[i + 8 * j] && (board[i + 8 * j].color == board.color || this.team())) {
			let ms = this.endFilter(board,1,i,j,this.moveFilter(board,i,j,this.getMoves(board,i,j)));
			ms.forEach((mv => {
				let b = this.doMove(board,i,j,mv);
				num += 1;
				this._propose(b,i0,j0,mv0,mv1,n - 1,zs2);
			}).bind(this));
		}
		
		let iv = window.setInterval((() => {
			if(zs2.length == num) {
				window.clearInterval(iv);
				let m0 = null;
				let xs = {};
				zs2.forEach(((m) => {
					if(m0 == null) {
						m0 = m[0];
						xs = {};
					}
					if(this.team()) {
						if(m[0] > m0) {
							m0 = m[0];
							xs = {};
						}

					} else {
						if(board.color == "white" && m[0] > m0) {
							m0 = m[0];
							xs = {};
						} else if(board.color == "black" && m[0] < m0) {

							m0 = m[0];
							xs = {};
						}
					}
					if(m[0] == m0) {
						xs[[i0,j0,mv0,mv1] + ""] = 1;
					}
				}).bind(this));
				zs.push([m0,xs]);
			}
		}).bind(this),100);
	}
	
	propose(n) {
		this.thinking = true;
		
		this.workers = [];
		for(let i = 0;i < 32;++i) this.workers.push(new Worker("ChessWorker.js"));
		this.workers.forEach(w => {
			w.done = true;
		});
		this.w = 0;

		let board = this.board;
		
		let num = 0;
		let zs = [];

		for(let i = 0;i < 8;++i)
		for(let j = 0;j < 8;++j)
		if(board[i + 8 * j] && (board[i + 8 * j].color == board.color || this.team())) {
			let ms = this.endFilter(board,1,i,j,this.moveFilter(board,i,j,this.getMoves(board,i,j)));
			ms.forEach((mv => {
				num += 1;
				let b = this.doMove(board,i,j,mv);
				this._propose(b,i,j,mv[0],mv[1],n,zs)
			}).bind(this));
		}
		this.progress = 1;
		this.render();		

		let iv = window.setInterval((() => {
			if(zs.length == num) {
				this.progress = 100.0 * zs.length / (num + 1);
				this.render();

				window.clearInterval(iv);
				
				let m0 = null;
				let xs = {};

				zs.forEach(((m) => {
					if(m0 == null) {
						m0 = m[0];
						xs = {};
					}
					if(this.team()) {
						if(m[0] > m0) {
							m0 = m[0];
							xs = {};
						}

					} else {
						if(board.color == "white" && m[0] > m0) {
							m0 = m[0];
							xs = {};
						} else if(board.color == "black" && m[0] < m0) {

							m0 = m[0];
							xs = {};
						}
					}
					if(m[0] == m0) {
						Object.keys(m[1]).forEach(y => {
							xs[y] = 1;
						});
					}
				}).bind(this));
				

				this.pp = [m0,xs];
				this.pp2 = {};
				Object.keys(this.pp[1]).forEach(p => {
					this.pp2[p.split(",").slice(0,2).join(",")] = 1;
				});

				console.log(this.pp);

				console.log(this.board);

				this.render();
				if(this.measureWhite)
					if(document.getElementById(this.cvID + "0"))
						document.getElementById(this.cvID + "0").innerHTML = "turn(" + this.board.color + ")  <br> points(" + this.measureWhite(this.board) + ")";
				this.progress = null;
				this.render();
				this.thinking = false;
			} else {
				this.progress = 100.0 * zs.length / (num + 1);
				this.render();
			}
		}).bind(this),100);
	}
}



class ClassicChess extends Chess {
	name() {
		return "Klassiches Schach";
	}
	
	measureWhite(board) {
		let cnt = 0;

		for(let i = 0;i < 8;++i)
		for(let j = 0;j < 8;++j)
		if(board[i + 8 * j] && board[i + 8 * j].color == "white")
			cnt += this.ps[board[i + 8 * j].shape];
		else
		if(board[i + 8 * j] && board[i + 8 * j].color == "black")
			cnt -= this.ps[board[i + 8 * j].shape];
		
		return cnt;
	}
	
	init() {
		this.board = [];
		for(var i = 0;i < 64;++i) this.board.push(null);
		
		this.board[0] = this.loadPiece("rook","rdd129","black");
		this.board[1] = this.loadPiece("knight","ndd129","black");
		this.board[2] = this.loadPiece("bishop","bdd129","black");
		this.board[3] = this.loadPiece("queen","qdd129","black");
		this.board[4] = this.loadPiece("king","kdd129","black");
		this.board[5] = this.loadPiece("bishop","bdd129","black");
		this.board[6] = this.loadPiece("knight","ndd129","black");
		this.board[7] = this.loadPiece("rook","rdd129","black");
		
		for(var i = 0;i < 8;++i) this.board[8 + i] = this.loadPiece("pawn","pdd129","black");
		
		this.board[56] = this.loadPiece("rook","rld129","white");
		this.board[57] = this.loadPiece("knight","nld129","white");
		this.board[58] = this.loadPiece("bishop","bld129","white");
		this.board[59] = this.loadPiece("queen","qld129","white");
		this.board[60] = this.loadPiece("king","kld129","white");
		this.board[61] = this.loadPiece("bishop","bld129","white");
		this.board[62] = this.loadPiece("knight","nld129","white");
		this.board[63] = this.loadPiece("rook","rld129","white");
		
		for(var i = 0;i < 8;++i) this.board[48 + i] = this.loadPiece("pawn","pld129","white");
	}

	moveFilter(board,i,j,xs) {
		var ys = [];
		var piece = board[i + 8 * j];
		
		for(var k = 0;k < xs.length;++k) {
			var x = xs[k];
			
			if(x[0] >= 0 && x[0] <= 7 && x[1] >= 0 && x[1] <= 7) {			
				switch(piece.shape) {
					case "pawn":
						if(
							(piece.color == "white" && j > x[1]) ||
							(piece.color == "black" && j < x[1])
						) {
							if(Math.abs(x[0] - i) == 1) {
								if(board[x[0] + 8 * x[1]] && board[x[0] + 8 * x[1]].color != piece.color)
									ys.push(x);
							} else if(!board[x[0] + 8 * x[1]]){
								if(piece.color == "black" && j == 1)
									ys.push(x);
								else if(piece.color == "white" && j == 6)
									ys.push(x);
								else if(Math.abs(x[1] - j) == 1)
									ys.push(x);
							}
						}
						break;
					case "rook":
					case "bishop":
					case "queen":
						var found = false;
						var l;
						for(l = 1;l < Math.max(Math.abs(i - x[0]),Math.abs(j - x[1]));++l)
							if(
								board[Math.sign(x[0] - i) * l + i + 8 * (l * Math.sign(x[1] - j) + j)]
								
							) {
								found = true;
								break;
							}
								
						if(!found)
							if(board[x[0] + 8 * x[1]] && board[x[0] + 8 * x[1]].color != piece.color)
								ys.push(x);
							else if(!board[x[0] + 8 * x[1]])
								ys.push(x);
						break;
					case "knight":
					case "king":
						if(board[x[0] + 8 * x[1]] && board[x[0] + 8 * x[1]].color != piece.color)
							ys.push(x);
						else if(!board[x[0] + 8 * x[1]])
							ys.push(x);
						break;
					default:
						ys.push(x);
				}
			}
			
		}
		
		var zs = [];
		for(var y = 0;y < ys.length;++y) { 
			var found = false;
			
			for(var z = 0;z < zs.length;++z) 
				if((ys[y][0] == zs[z][0]) && (ys[y][1] == zs[z][1]))
					found = true;
			
			if(!found) zs.push(ys[y]);
		}
		
		return zs;
	}
	
	endFilter(board,n,i,j,xs) {
		if(n == 0) return xs;

		var ys = [];
		var zs = board[i + 8 * j].color == "white" ? ["king"] : ["king"];
		
		for(var k = 0;k < xs.length;++k) {
			var found = false;
			var x = xs[k];
			var brd = this.doMove(board,i,j,x);
			for(var h = 0;h < zs.length;++h) {
				for(var a = 0;a < 8;++a) {
					for(var b = 0;b < 8;++b) {
						if(brd[a + 8 * b] && (brd[a + 8 * b].shape == zs[h]) && (board[i + 8 * j].color == brd[a + 8 * b].color)) {
							if(this.isUnderThreat(brd,n,a,b)) found = true;
						}
					}
				}
			}
			if(!found) ys.push(x);
		}
		
		zs = [];
		for(var y = 0;y < ys.length;++y) { 
			var found = false;
			
			for(var z = 0;z < zs.length;++z) 
				if((ys[y][0] == zs[z][0]) && (ys[y][1] == zs[z][1]))
					found = true;
			
			if(!found) zs.push(ys[y]);
		}
		
		return zs;
	}
	
	_doMove(board,i,j,m) {	
		var brd = [];
		for(var k = 0;k < 64;++k) brd.push(board[k]);
		
		if(
			m[0] == 6 && m[1] == 0
			&& i == 4 && j == 0 && brd[4] && brd[4].color == "black" && brd[4].shape == "king" 
			&& !brd[5] && !brd[6] && brd[7] && brd[7].color == "black" && brd[7].shape == "rook"
		) {
			brd[5] = brd[7];
			brd[7] = null;
		}
		
		if(
			m[0] == 2 && m[1] == 0
			&& i == 4 && j == 0 && brd[4] && brd[4].color == "black" && brd[4].shape == "king" 
			&& !brd[3] && !brd[2] && !brd[1] && brd[0] && brd[0].color == "black" && brd[0].shape == "rook"
		) {
			brd[3] = brd[0];
			brd[0] = null;
		}
		
		
		if(
			m[0] == 6 && m[1] == 7
			&& i == 4 && j == 7 && brd[60] && brd[60].color == "white" && brd[60].shape == "king" 
			&& !brd[61] && !brd[62] && brd[63] && brd[63].color == "white" && brd[63].shape == "rook"
		) {
			brd[61] = brd[63];
			brd[63] = null;
		}
		
		if(
			m[0] == 2 && m[1] == 7
			&& i == 4 && j == 7 && brd[60] && brd[60].color == "white" && brd[60].shape == "king" 
			&& !brd[59] && !brd[58] && !brd[57] && brd[56] && brd[56].color == "white" && brd[56].shape == "rook"
		) {
			brd[59] = brd[56];
			brd[56] = null;
		}

		brd[m[0] + 8 * m[1]] = brd[i + 8 * j];
		brd[i + 8 * j] = null;

		return brd;
	}
}



class ZauberStrikeChess extends ClassicChess {
	name() {
		return "Zauber Strike Schach";
	}

	measureWhite(board) {
		let cnt = 0;

		for(let i = 0;i < 8;++i)
		for(let j = 0;j < 8;++j)
		if(board[i + 8 * j])
			cnt += this.ps[board[i + 8 * j].shape];
		
		return 32 - cnt;
	}

	init() {
		this.ps["queen"] = 1;
		this.ps["rook"] = 1;
		this.ps["knight"] = 1;
		this.ps["bishop"] = 1;
		this.ps["pawn"] = 1;
		this.ps["king"] = 1;
		
		this.board = [];
		for(var i = 0;i < 64;++i) this.board.push(null);
		
		this.board[0] = this.loadPiece("rook","rdd129","black");
		this.board[1] = this.loadPiece("knight","ndd129","black");
		this.board[2] = this.loadPiece("bishop","bdd129","black");
		this.board[3] = this.loadPiece("queen","qdd129","black");
		this.board[4] = this.loadPiece("king","kdd129","black");
		this.board[5] = this.loadPiece("bishop","bdd129","black");
		this.board[6] = this.loadPiece("knight","ndd129","black");
		this.board[7] = this.loadPiece("rook","rdd129","black");
		
		for(var i = 0;i < 8;++i) this.board[8 + i] = this.loadPiece("pawn","pdd129","black");
		
		this.board[56] = this.loadPiece("rook","rld129","white");
		this.board[57] = this.loadPiece("knight","nld129","white");
		this.board[58] = this.loadPiece("bishop","bld129","white");
		this.board[59] = this.loadPiece("queen","qld129","white");
		this.board[60] = this.loadPiece("king","kld129","white");
		this.board[61] = this.loadPiece("bishop","bld129","white");
		this.board[62] = this.loadPiece("knight","nld129","white");
		this.board[63] = this.loadPiece("rook","rld129","white");
		
		for(var i = 0;i < 8;++i) this.board[48 + i] = this.loadPiece("pawn","pld129","white");
		this.board = this.shuffle(this.board);
	}

	moveFilter(board,i,j,xs) {
		var ys = [];
		var piece = board[i + 8 * j];
		
		for(var k = 0;k < xs.length;++k) {
			var x = xs[k];
			
			if(x[0] >= 0 && x[0] <= 7 && x[1] >= 0 && x[1] <= 7) {			
				switch(piece.shape) {
					case "pawn":
						break;
					default:
						if(board[x[0] + 8 * x[1]] && board[x[0] + 8 * x[1]].color != piece.color)
							ys.push(x);
						break;
				}
			}
			
		}
		
		var zs = [];
		for(var y = 0;y < ys.length;++y) { 
			var found = false;
			
			for(var z = 0;z < zs.length;++z) 
				if((ys[y][0] == zs[z][0]) && (ys[y][1] == zs[z][1]))
					found = true;
			
			if(!found) zs.push(ys[y]);
		}
		
		return zs;
	}

	endFilter(board,n,i,j,xs) {
		return xs;
	}

	_doMove(board,i,j,m) {	
		var brd = [];
		for(var k = 0;k < 64;++k) brd.push(board[k]);
		brd[m[0] + 8 * m[1]] = brd[i + 8 * j];
		brd[i + 8 * j] = null;

		return brd;
	}
}











class GardenEdenChess extends Chess {
	name() {
		return "Garten Eden Schach";
	}
	
	measureWhite(board) {
		let col = board[0].color
		let cnt = 0;
		
		for(let i = 0;i < 8;++i) {
			if(board[i].color == col) cnt += 1;
			col = col == "white" ? "black" : "white";
		}
		for(let i = 1;i < 7;++i) {
			if(board[7 + 8 * i].color == col) cnt += 1;
			col = col == "white" ? "black" : "white";
		}
		for(let i = 0;i < 8;++i) {
			if(board[63 - i].color == col) cnt += 1;
			col = col == "white" ? "black" : "white";
		}
		for(let i = 1;i < 7;++i) {
			if(board[0 + 8 * (7 - i)].color == col) cnt += 1;
			col = col == "white" ? "black" : "white";
		}
		
		for(let i = 3;i <= 4;++i) {
			for(let j = 3;j <= 4;++j) {
				if(board[i + 8 * j].shape == "king" || board[i + 8 * j].shape == "queen")
					cnt += 1;
			}
		}
		
		return Math.round(100.0 * cnt / 32.0);
	}


	init() {
		this.board = [];
		for(var i = 0;i < 64;++i) this.board.push(null);
		
		var brd = [];
		
		brd.push(this.loadPiece("rook","rdd129","black"));
		brd.push(this.loadPiece("knight","ndd129","black"));
		brd.push(this.loadPiece("bishop","bdd129","black"));
		brd.push(this.loadPiece("queen","qdd129","black"));
		brd.push(this.loadPiece("king","kdd129","black"));
		brd.push(this.loadPiece("bishop","bdd129","black"));
		brd.push(this.loadPiece("knight","ndd129","black"));
		brd.push(this.loadPiece("rook","rdd129","black"));
		for(var i = 0;i < 8;++i) brd.push(this.loadPiece("pawn","pdd129","black"));

		brd.push(this.loadPiece("rook","rld129","white"));
		brd.push(this.loadPiece("knight","nld129","white"));
		brd.push(this.loadPiece("bishop","bld129","white"));
		brd.push(this.loadPiece("queen","qld129","white"));
		brd.push(this.loadPiece("king","kld129","white"));
		brd.push(this.loadPiece("bishop","bld129","white"));
		brd.push(this.loadPiece("knight","nld129","white"));
		brd.push(this.loadPiece("rook","rld129","white"));
		for(var i = 0;i < 8;++i) brd.push(this.loadPiece("pawn","pld129","white"));
		
		//brd = shuffle(brd);
		
		for(var i = 0;i < 8;++i) this.board[i] = brd.pop();
		for(var i = 0;i < 8;++i) this.board[i + 56] = brd.pop();
		for(var i = 0;i < 6;++i) this.board[8 + 8 * i] = brd.pop();
		for(var i = 0;i < 6;++i) this.board[15 + 8 * i] = brd.pop();
		
		this.board[27] = brd.pop();
		this.board[28] = brd.pop();
		this.board[35] = brd.pop();
		this.board[36] = brd.pop();
	}
	

	moveFilter(board,i,j,xs) {
		var ys = [];
		var piece = board[i + 8 * j];
		
		for(var k = 0;k < xs.length;++k) {
			var x = xs[k];
			
			if(x[0] >= 0 && x[0] <= 7 && x[1] >= 0 && x[1] <= 7) {			
				switch(piece.shape) {
					case "pawn":
						break;
					case "rook":
					case "bishop":
					case "queen":
						var found = false;
						var l;
						for(l = 1;l < Math.max(Math.abs(i - x[0]),Math.abs(j - x[1]));++l)
							if(
								board[Math.sign(x[0] - i) * l + i + 8 * (l * Math.sign(x[1] - j) + j)]
								
							) {
								found = true;
								break;
							}
								
						if(!found)
							if(board[x[0] + 8 * x[1]])
								ys.push(x);
						break;
					case "knight":
					case "king":
						if(board[x[0] + 8 * x[1]])
							ys.push(x);
						break;
					default:
						ys.push(x);
				}
			}
			
		}
		
		var zs = [];
		for(var y = 0;y < ys.length;++y) { 
			var found = false;
			
			for(var z = 0;z < zs.length;++z) 
				if((ys[y][0] == zs[z][0]) && (ys[y][1] == zs[z][1]))
					found = true;
			
			if(!found) zs.push(ys[y]);
		}
		
		return zs;
	}

	endFilter(board,n,i,j,xs) {
		return xs;
	}

	_doMove(board,i,j,m) {
		var brd = [];
		for(var k = 0;k < 64;++k) brd.push(board[k]);
		
		var swp = brd[m[0] + 8 * m[1]];
		brd[m[0] + 8 * m[1]] = brd[i + 8 * j];
		brd[i + 8 * j] = swp;

		return brd;
	}
}














class HalmaChess extends GardenEdenChess {
	name() {
		return "Halma Schach";
	}
	
	measureWhite(board) {
		let cnt = 0;
		
		if(board[3 + 4 * 8] && board[3 + 4 * 8].color == "white") cnt += 1;			
		for(let i = 0;i <= 3;++i)
			for(let j = 0;j <= i;++j)
				if(board[4 * 8 + j + i * 8] &&
					board[4 * 8 + j + 8 * i].color == "white" &&
						board[4 * 8 + j + 8 * i].shape != "pawn"
				) cnt += 1;
			
		for(let i = 0;i < 5;++i)
			if(board[3 * 8 + 9 * i] &&
				board[3 * 8 + 9 * i].color == "white" &&
					board[3 * 8 + 9 * i].shape == "pawn"
			) cnt += 1;
	
		if(board[4 + 3 * 8] && board[4 + 3 * 8].color == "black") cnt -= 1;			
		for(let i = 0;i <= 3;++i)
			for(let j = 0;j <= i;++j)
				if(board[4 + i + j * 8] &&
					board[4 + i + 8 * j].color == "black" &&
						board[4 + i + 8 * j].shape != "pawn"
				) cnt -= 1;
			
		for(let i = 0;i < 5;++i)
			if(board[3 + 9 * i] &&
				board[3 + 9 * i].color == "black" &&
					board[3 + 9 * i].shape == "pawn"
			) cnt -= 1;

		return cnt;
	}

	init() {
		this.board = [];
		for(var i = 0;i < 64;++i) this.board.push(null);
		
		var brd = [];
		var brd2 = [];
		
		brd.push(this.loadPiece("rook","rdd129","black"));
		brd.push(this.loadPiece("knight","ndd129","black"));
		brd.push(this.loadPiece("bishop","bdd129","black"));
		brd.push(this.loadPiece("queen","qdd129","black"));
		brd.push(this.loadPiece("king","kdd129","black"));
		brd.push(this.loadPiece("bishop","bdd129","black"));
		brd.push(this.loadPiece("knight","ndd129","black"));
		brd.push(this.loadPiece("rook","rdd129","black"));
		for(var i = 0;i < 8;++i) brd.push(this.loadPiece("pawn","pdd129","black"));

		brd2.push(this.loadPiece("rook","rld129","white"));
		brd2.push(this.loadPiece("knight","nld129","white"));
		brd2.push(this.loadPiece("bishop","bld129","white"));
		brd2.push(this.loadPiece("queen","qld129","white"));
		brd2.push(this.loadPiece("king","kld129","white"));
		brd2.push(this.loadPiece("bishop","bld129","white"));
		brd2.push(this.loadPiece("knight","nld129","white"));
		brd2.push(this.loadPiece("rook","rld129","white"));
		for(var i = 0;i < 8;++i) brd2.push(this.loadPiece("pawn","pld129","white"));
		
		brd = shuffle(brd);
		brd2 = shuffle(brd2);
		
		for(var i = 0;i <= 4;++i) for(var j = 0;j <= i;++j) this.board[(4 - i) + 8 * (7 - j)] = brd.pop();
		for(var i = 0;i <= 4;++i) for(var j = 0;j <= i;++j) this.board[(i + 3) + 8 * j] = brd2.pop();
		this.board[3 + 8 * 4] = brd.pop();
		this.board[4 + 8 * 3] = brd2.pop();
	}
	

	moveFilter(board,i,j,xs) {
		var ys = [];
		var piece = board[i + 8 * j];
		
		for(var k = 0;k < xs.length;++k) {
			var x = xs[k];
			
			if(x[0] >= 0 && x[0] <= 7 && x[1] >= 0 && x[1] <= 7) {			
				switch(piece.shape) {
					case "pawn":
						ys.push(x);
						break;
					case "rook":
					case "bishop":
					case "queen":
						var found = false;
						var l;
						for(l = 1;l < Math.max(Math.abs(i - x[0]),Math.abs(j - x[1]));++l)
							if(
								board[Math.sign(x[0] - i) * l + i + 8 * (l * Math.sign(x[1] - j) + j)]
								
							) {
								found = true;
								break;
							}
								
						if(!found)
								ys.push(x);
						break;
					case "knight":
					case "king":
							ys.push(x);
						break;
					default:
						ys.push(x);
				}
			}
			
		}
		
		var zs = [];
		for(var y = 0;y < ys.length;++y) { 
			var found = false;
			
			for(var z = 0;z < zs.length;++z) 
				if((ys[y][0] == zs[z][0]) && (ys[y][1] == zs[z][1]))
					found = true;
			
			if(!found) zs.push(ys[y]);
		}
		
		return zs;
	}

}





class GalagaChess extends HalmaChess {
	name() {
		return "Galaga Schach";
	}
	
	init() {
		this.board = [];
		for(var i = 0;i < 64;++i) this.board.push(null);
		
		var brd = [];
		var brd2 = [];
		
		for(var i = 0;i < 8;++i) this.board[i + 8 * i] = this.loadPiece("pawn","pdd129","black");
		for(var i = 0;i < 8;++i) this.board[(7 - i) + 8 * i] = this.loadPiece("pawn","pld129","white");
		
		brd.push(this.loadPiece("queen","qdd129","black"));
		brd.push(this.loadPiece("king","kdd129","black"));
		brd.push(this.loadPiece("queen","qld129","white"));
		brd.push(this.loadPiece("king","kld129","white"));
		for(let i = 0;i < 8;++i) brd.push(null);

		brd2.push(this.loadPiece("rook","rdd129","black"));
		brd2.push(this.loadPiece("knight","ndd129","black"));
		brd2.push(this.loadPiece("bishop","bdd129","black"));
		brd2.push(this.loadPiece("bishop","bdd129","black"));
		brd2.push(this.loadPiece("knight","ndd129","black"));
		brd2.push(this.loadPiece("rook","rdd129","black"));
		brd2.push(this.loadPiece("rook","rld129","white"));
		brd2.push(this.loadPiece("knight","nld129","white"));
		brd2.push(this.loadPiece("bishop","bld129","white"));
		brd2.push(this.loadPiece("bishop","bld129","white"));
		brd2.push(this.loadPiece("knight","nld129","white"));
		brd2.push(this.loadPiece("rook","rld129","white"));
		
		brd = shuffle(brd);
		brd2 = shuffle(brd2);
		
		for(var i = 0;i <= 2;++i) {
			for(var j = 0;j <= i;++j) this.board[(i + 1) + 8 * (7 - j)] = brd2.pop();
			for(var j = 0;j <= i;++j) this.board[(6 - i) + 8 * (7 - j)] = brd2.pop();
		}
		for(var i = 0;i <= 2;++i) {
			for(var j = 0;j <= i;++j) this.board[(i + 1) + 8 * j] = brd.pop();
			for(var j = 0;j <= i;++j) this.board[(6 - i) + 8 * j] = brd.pop();
		}
	}
	
	moveFilter(board,i,j,xs) {
		let ys = [];
		var piece = board[i + 8 * j];
		
		xs.forEach(x => {
			if(!(i > j && i < 7 - j))
			if(piece.shape != "pawn")
			if(x[0] >= 0 && x[0] <= 7 && x[1] >= 0 && x[1] <= 7) {
				if(x[1] < 4) {
					if(x[0] == x[1])
						;
					else if(x[0] == 7 - x[1])
						;
					else if(x[0] > x[1] && x[0] < 7 - x[1]) {
						if(board[x[0] + 8 * x[1]])
							ys.push(x);
					} else
						ys.push(x);
				} else {
					ys.push(x);
				}
			}
		});
		
		let zs = [];
		ys.forEach(x => {
			switch(piece.shape) {
				case "pawn":
					break;
				case "rook":
				case "bishop":
				case "queen":
					var found = false;
					var l;
					for(l = 1;l < Math.max(Math.abs(i - x[0]),Math.abs(j - x[1]));++l) {
						let p = board[Math.sign(x[0] - i) * l + i + 8 * (l * Math.sign(x[1] - j) + j)];
						if(l * Math.sign(x[1] - j) + j >= 4)
						if(p && p.shape == "pawn") {
							found = true;
							break;
						}
					}
							
					if(!found)
							zs.push(x);
					break;
				case "knight":
				case "king":
						zs.push(x);
					break;
				default:
					zs.push(x);
			}
		});
		
		let cnt = 0;
		
		for(let i = 0;i < 8;++i)
		for(let j = 0;j < 8;++j)
			if(board[i + 8 * j] && board[i + 8 * j].shape != "pawn" && board[i + 8 * j].color == board.color)
				if(j < 4) {
					if(i > j && i < 7 - j) {
					} else
						cnt += 1;
				} else {
					cnt += 1;
				}
		
		if(cnt < 6) {
			let bs = [];
			zs.forEach(x => {
				if(x[0] > x[1] && x[0] < 7 - x[1])
					if(board[x[0] + 8 * x[1]])
						bs.push(x);
			});
			
			return bs;
		} else {
			if(!this.recursion) {
				let c = 0;

				for(let i = 0;i < 8;++i)
				for(let j = 0;j < 8;++j)
				if(board[i + 8 * j] && board[i + 8 * j].shape != "pawn" && board[i + 8 * j].color == board.color) {
					this.recursion = true;
					let ms = this.moveFilter(board,i,j,this.getMoves(board,i,j));
					this.recursion = false;
					ms.forEach(m => {
						if(m[0] > m[1] && m[0] < 7 - m[1]) c += 1;
					});
				}
				
				if(c < 6) {
					let bs = [];
					zs.forEach(x => {
						if(x[0] > x[1] && x[0] < 7 - x[1]) 
							if(board[x[0] + 8 * x[1]])
								return;
						bs.push(x);
					});
					
					return bs;
				}
			}
		}
		
		return zs;
	}
	
	//next board with the former color
	next(board) {
		let c = board.color == "white" ? "black" : "white";

		for(let i = 0;i < 8;++i)
		for(let j = 0;j < 8;++j)
			if(board[i + 8 * j] && board[i + 8 * j].shape != "pawn" && board[i + 8 * j].color == c) {
				let ms = this.moveFilter(board,i,j,this.getMoves(board,i,j));
				if(ms.length > 0) return c;
			}
		
		return board.color;
	}

	_doMove(board,i,j,m) {	
		var brd = [];
		for(var k = 0;k < 64;++k) brd.push(board[k]);
		let swp = brd[m[0] + 8 * m[1]];
		brd[m[0] + 8 * m[1]] = brd[i + 8 * j];
		if(m[1] < 4) {
			if(m[0] > m[1] && m[0] < 7 - m[1]) {
				brd[i + 8 * j] = null;
			} else {
				brd[i + 8 * j] = swp;
			}
		} else {
			brd[i + 8 * j] = swp;
		}

		return brd;
	}
	
	measureWhite(board) {
		let cnt = 0;
		
		
		for(let i = 0;i < 8;++i)
		for(let j = 0;j < 8;++j)
			if(board[i + 8 * j] && board[i + 8 * j].shape != "pawn")
				if(j < 4) {
					if(i > j && i < 7 - j) {
					} else
						cnt += 2 * (board[i + 8 * j].color == "white" ? 1 : -1);
				} else {
					cnt += 2 * (board[i + 8 * j].color == "white" ? 1 : -1);
				}
		cnt = -cnt;
		
		for(let i = 0;i < 8;++i)
		for(let j = 0;j < 8;++j)
			if(board[i + 8 * j] && board[i + 8 * j].shape != "pawn")
				if(i > j && i < 7 - j) {
				} else {
					let ms = this.moveFilter(board,i,j,this.getMoves(board,i,j));
					for(let k = 0;k < ms.length;++k) {
						let m = ms[k];
						if(m[0] > m[1] && m[0] < 7 - m[1]) {
							cnt +=  1 * (board[i + 8 * j].color == "white" ? 1 : -1);
							break;
						}
					}
				}
		return cnt;
	}
}





class MonopolyChess extends HalmaChess {
	name() {
		return "Monopoly Schach";
	}
	
	isRing1(board,i,j) {
		if(
			(i == 3 && j >= 3 && j <= 4) ||
			(i == 4 && j >= 3 && j <= 4) ||
			(j == 3 && i >= 3 && i <= 4) ||
			(j == 4 && i >= 3 && i <= 4) 
		)
		if(board[i + 8 * j]) {
			if(board[i + 8 * j].shape == "king") {
				return true;
			}
			if(board[i + 8 * j].shape == "queen") {
				return true;
			}
		}
		return false;
	}

	isRing2(board,i,j) {
		if(
			(i == 2 && j >= 2 && j <= 5) ||
			(i == 5 && j >= 2 && j <= 5) ||
			(j == 3 && i >= 2 && i <= 5) ||
			(j == 5 && i >= 2 && i <= 5) 
		)
		if(board[i + 8 * j]) {
			if(board[i + 8 * j].shape == "knight") {
				return true;
			}
			if(board[i + 8 * j].shape == "bishop") {
				return true;
			}
			if(board[i + 8 * j].shape == "rook") {
				return true;
			}
		}
		return false;
	}
	
	isRing3(board,i,j) {
		if(
			(i == 1 && j >= 1 && j <= 6) ||
			(i == 6 && j >= 1 && j <= 6) ||
			(j == 1 && i >= 1 && i <= 6) ||
			(j == 6 && i >= 1 && i <= 6) 
		)
		if(board[i + 8 * j]) {
			if(board[i + 8 * j].shape == "pawn") {
				return true;
			}
		}
		return false;
	}
	
	measureWhite(board) {
		let cnt = 0;
		
		for(let i = 0;i < 8;++i) {
			for(let j = 0;j < 8;++j) {
				if(this.isRing1(board,i,j)) cnt += 1;
				else if(this.isRing2(board,i,j)) cnt += 1;
				else if(this.isRing3(board,i,j)) cnt += 1;
			}
		}

		return Math.round(1000.0 * cnt / (board.counter + 1));
	}
	
	init() {
		this.board = [];
		for(var i = 0;i < 64;++i) this.board.push(null);
		
		var brd = [];
		
		brd.push(this.loadPiece("rook","rdd129","black"));
		brd.push(this.loadPiece("knight","ndd129","black"));
		brd.push(this.loadPiece("bishop","bdd129","black"));
		brd.push(this.loadPiece("queen","qdd129","black"));
		brd.push(this.loadPiece("king","kdd129","black"));
		brd.push(this.loadPiece("bishop","bdd129","black"));
		brd.push(this.loadPiece("knight","ndd129","black"));
		brd.push(this.loadPiece("rook","rdd129","black"));
		for(var i = 0;i < 8;++i) brd.push(this.loadPiece("pawn","pdd129","black"));

		brd.push(this.loadPiece("rook","rld129","white"));
		brd.push(this.loadPiece("knight","nld129","white"));
		brd.push(this.loadPiece("bishop","bld129","white"));
		brd.push(this.loadPiece("queen","qld129","white"));
		brd.push(this.loadPiece("king","kld129","white"));
		brd.push(this.loadPiece("bishop","bld129","white"));
		brd.push(this.loadPiece("knight","nld129","white"));
		brd.push(this.loadPiece("rook","rld129","white"));
		for(var i = 0;i < 8;++i) brd.push(this.loadPiece("pawn","pld129","white"));
		
		brd = shuffle(brd);
		
		for(var i = 0;i < 8;++i) this.board[i] = brd.pop();
		for(var i = 0;i < 8;++i) this.board[i + 56] = brd.pop();
		for(var i = 0;i < 6;++i) this.board[8 + 8 * i] = brd.pop();
		for(var i = 0;i < 6;++i) this.board[15 + 8 * i] = brd.pop();
		
		this.board[9] = brd.pop();
		this.board[14] = brd.pop();
		this.board[49] = brd.pop();
		this.board[54] = brd.pop();
	}
	

	moveFilter(board,i,j,xs) {
		let ys = [];
		var piece = board[i + 8 * j];
		
		xs.forEach(x => {
			if(x[0] >= 0 && x[0] <= 7 && x[1] >= 0 && x[1] <= 7) {
				if(piece.shape != "pawn")
					ys.push(x);
			}
		});
				
		return ys;
	}

	endFilter(board,n,i,j,xs) {
		return xs;
	}

	_doMove(board,i,j,m) {
		var brd = [];
		for(var k = 0;k < 64;++k) brd.push(board[k]);
		
		var swp = brd[m[0] + 8 * m[1]];
		brd[m[0] + 8 * m[1]] = brd[i + 8 * j];
		brd[i + 8 * j] = swp;

		return brd;
	}
}



function md5cycle(x, k) {
var a = x[0], b = x[1], c = x[2], d = x[3];

a = ff(a, b, c, d, k[0], 7, -680876936);
d = ff(d, a, b, c, k[1], 12, -389564586);
c = ff(c, d, a, b, k[2], 17,  606105819);
b = ff(b, c, d, a, k[3], 22, -1044525330);
a = ff(a, b, c, d, k[4], 7, -176418897);
d = ff(d, a, b, c, k[5], 12,  1200080426);
c = ff(c, d, a, b, k[6], 17, -1473231341);
b = ff(b, c, d, a, k[7], 22, -45705983);
a = ff(a, b, c, d, k[8], 7,  1770035416);
d = ff(d, a, b, c, k[9], 12, -1958414417);
c = ff(c, d, a, b, k[10], 17, -42063);
b = ff(b, c, d, a, k[11], 22, -1990404162);
a = ff(a, b, c, d, k[12], 7,  1804603682);
d = ff(d, a, b, c, k[13], 12, -40341101);
c = ff(c, d, a, b, k[14], 17, -1502002290);
b = ff(b, c, d, a, k[15], 22,  1236535329);

a = gg(a, b, c, d, k[1], 5, -165796510);
d = gg(d, a, b, c, k[6], 9, -1069501632);
c = gg(c, d, a, b, k[11], 14,  643717713);
b = gg(b, c, d, a, k[0], 20, -373897302);
a = gg(a, b, c, d, k[5], 5, -701558691);
d = gg(d, a, b, c, k[10], 9,  38016083);
c = gg(c, d, a, b, k[15], 14, -660478335);
b = gg(b, c, d, a, k[4], 20, -405537848);
a = gg(a, b, c, d, k[9], 5,  568446438);
d = gg(d, a, b, c, k[14], 9, -1019803690);
c = gg(c, d, a, b, k[3], 14, -187363961);
b = gg(b, c, d, a, k[8], 20,  1163531501);
a = gg(a, b, c, d, k[13], 5, -1444681467);
d = gg(d, a, b, c, k[2], 9, -51403784);
c = gg(c, d, a, b, k[7], 14,  1735328473);
b = gg(b, c, d, a, k[12], 20, -1926607734);

a = hh(a, b, c, d, k[5], 4, -378558);
d = hh(d, a, b, c, k[8], 11, -2022574463);
c = hh(c, d, a, b, k[11], 16,  1839030562);
b = hh(b, c, d, a, k[14], 23, -35309556);
a = hh(a, b, c, d, k[1], 4, -1530992060);
d = hh(d, a, b, c, k[4], 11,  1272893353);
c = hh(c, d, a, b, k[7], 16, -155497632);
b = hh(b, c, d, a, k[10], 23, -1094730640);
a = hh(a, b, c, d, k[13], 4,  681279174);
d = hh(d, a, b, c, k[0], 11, -358537222);
c = hh(c, d, a, b, k[3], 16, -722521979);
b = hh(b, c, d, a, k[6], 23,  76029189);
a = hh(a, b, c, d, k[9], 4, -640364487);
d = hh(d, a, b, c, k[12], 11, -421815835);
c = hh(c, d, a, b, k[15], 16,  530742520);
b = hh(b, c, d, a, k[2], 23, -995338651);

a = ii(a, b, c, d, k[0], 6, -198630844);
d = ii(d, a, b, c, k[7], 10,  1126891415);
c = ii(c, d, a, b, k[14], 15, -1416354905);
b = ii(b, c, d, a, k[5], 21, -57434055);
a = ii(a, b, c, d, k[12], 6,  1700485571);
d = ii(d, a, b, c, k[3], 10, -1894986606);
c = ii(c, d, a, b, k[10], 15, -1051523);
b = ii(b, c, d, a, k[1], 21, -2054922799);
a = ii(a, b, c, d, k[8], 6,  1873313359);
d = ii(d, a, b, c, k[15], 10, -30611744);
c = ii(c, d, a, b, k[6], 15, -1560198380);
b = ii(b, c, d, a, k[13], 21,  1309151649);
a = ii(a, b, c, d, k[4], 6, -145523070);
d = ii(d, a, b, c, k[11], 10, -1120210379);
c = ii(c, d, a, b, k[2], 15,  718787259);
b = ii(b, c, d, a, k[9], 21, -343485551);

x[0] = add32(a, x[0]);
x[1] = add32(b, x[1]);
x[2] = add32(c, x[2]);
x[3] = add32(d, x[3]);

}

function cmn(q, a, b, x, s, t) {
a = add32(add32(a, q), add32(x, t));
return add32((a << s) | (a >>> (32 - s)), b);
}

function ff(a, b, c, d, x, s, t) {
return cmn((b & c) | ((~b) & d), a, b, x, s, t);
}

function gg(a, b, c, d, x, s, t) {
return cmn((b & d) | (c & (~d)), a, b, x, s, t);
}

function hh(a, b, c, d, x, s, t) {
return cmn(b ^ c ^ d, a, b, x, s, t);
}

function ii(a, b, c, d, x, s, t) {
return cmn(c ^ (b | (~d)), a, b, x, s, t);
}

function md51(s) {
txt = '';
var n = s.length,
state = [1732584193, -271733879, -1732584194, 271733878], i;
for (i=64; i<=s.length; i+=64) {
md5cycle(state, md5blk(s.substring(i-64, i)));
}
s = s.substring(i-64);
var tail = [0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0];
for (i=0; i<s.length; i++)
tail[i>>2] |= s.charCodeAt(i) << ((i%4) << 3);
tail[i>>2] |= 0x80 << ((i%4) << 3);
if (i > 55) {
md5cycle(state, tail);
for (i=0; i<16; i++) tail[i] = 0;
}
tail[14] = n*8;
md5cycle(state, tail);
return state;
}

/* there needs to be support for Unicode here,
 * unless we pretend that we can redefine the MD-5
 * algorithm for multi-byte characters (perhaps
 * by adding every four 16-bit characters and
 * shortening the sum to 32 bits). Otherwise
 * I suggest performing MD-5 as if every character
 * was two bytes--e.g., 0040 0025 = @%--but then
 * how will an ordinary MD-5 sum be matched?
 * There is no way to standardize text to something
 * like UTF-8 before transformation; speed cost is
 * utterly prohibitive. The JavaScript standard
 * itself needs to look at this: it should start
 * providing access to strings as preformed UTF-8
 * 8-bit unsigned value arrays.
 */
function md5blk(s) { /* I figured global was faster.   */
var md5blks = [], i; /* Andy King said do it this way. */
for (i=0; i<64; i+=4) {
md5blks[i>>2] = s.charCodeAt(i)
+ (s.charCodeAt(i+1) << 8)
+ (s.charCodeAt(i+2) << 16)
+ (s.charCodeAt(i+3) << 24);
}
return md5blks;
}

var hex_chr = '0123456789abcdef'.split('');

function rhex(n)
{
var s='', j=0;
for(; j<4; j++)
s += hex_chr[(n >> (j * 8 + 4)) & 0x0F]
+ hex_chr[(n >> (j * 8)) & 0x0F];
return s;
}

function hex(x) {
for (var i=0; i<x.length; i++)
x[i] = rhex(x[i]);
return x.join('');
}

function md5(s) {
return hex(md51(s));
}

/* this function is much faster,
so if possible we use it. Some IEs
are the only ones I know of that
need the idiotic second function,
generated by an if clause.  */

function add32(a, b) {
return (a + b) & 0xFFFFFFFF;
}

if (md5('hello') != '5d41402abc4b2a76b9719d911017c592') {
function add32(x, y) {
var lsw = (x & 0xFFFF) + (y & 0xFFFF),
msw = (x >> 16) + (y >> 16) + (lsw >> 16);
return (msw << 16) | (lsw & 0xFFFF);
}
}



function xmur3(str) {
    for(var i = 0, h = 1779033703 ^ str.length; i < str.length; i++)
        h = Math.imul(h ^ str.charCodeAt(i), 3432918353),
        h = h << 13 | h >>> 19;
    return function() {
        h = Math.imul(h ^ h >>> 16, 2246822507);
        h = Math.imul(h ^ h >>> 13, 3266489909);
        return (h ^= h >>> 16) >>> 0;
    }
}

function mulberry32(a) {
    return function() {
      var t = a += 0x6D2B79F5;
      t = Math.imul(t ^ t >>> 15, t | 1);
      t ^= t + Math.imul(t ^ t >>> 7, t | 61);
      return ((t ^ t >>> 14) >>> 0) / 4294967296;
    }
}




function Random(str) {
	var seed = xmur3(str);
	this.rnd = mulberry32(seed());
}
	



function Tarot(str) {
	this.rand = new Random(str);
	this.installedVerb = "";
}

Tarot.prototype.rnd = function() {
	return Math.round(this.rand.rnd() * 10000);
}
	
Tarot.prototype.clichee = function() {
	var database = [
		["nördlich","südlich","geh","in","einer Woche"],
		["östlich","westlich","geh","in","einer Woche"],
		["ober","unter","flieg","mit","dem Strom"],
		["link","recht","dreh","mit","dem Uhrzeigersinn"],
		["klein","groß","wachs","mit","den Jahreszeiten"],
		["hell","dunkel","verblass","in","Sekunden"],
		["fest","porös","verringer","mit","der Zeit"],
		["rund","eckig","verform","nach","der Formel"],
		["weiblich","männlich","transformier","mit","der Therapie"],
		["gut","boese","konvertier","mit","der Sünde"],
		["alt","jung"],
		["alt","neu"],
		["alt","neu"],
		["heiß","kalt"],
	];
	
	var db = [
		["in","einer Woche"],
		["in","einer Woche"],
		["mit","dem Strom"],
		["mit","dem Uhrzeigersinn"],
		["mit","den Jahreszeiten"],
		["in","einem Jahr"],
		["in","Sekunden"],
		["mit","der Zeit"],
		["nach","der Formel"],
		["mit","der Therapie"],
		["mit","der Sünde"],
	];

	let a = database[this.rnd() % database.length];
	let b = db[this.rnd() % db.length];
	
	return [a[0],a[1],this.verb(),b[0],b[1]];
}
	
Tarot.prototype.verb = function() {
	var words = [
	
		"back",
		"fahr",
		"grab",
		"lad",
		"schaff",
		"schlag",
		"trag",
		"wachs",
		"wasch",
		"blas",
		"brat",
		"fall",
		"halt",
		"lass",
		"rat",
		"schlaf",
		"empfang",
		"empfang",
		"fang",
		"les",
		"seh",
		"befehl",
		"befohl",
		"empfehl",
		"empfohl",
		"stehl",
		"ess",
		"fress",
		"mess",
		"tret",
		"vergess",
		"vergess",
		"berg",
		"berst",
		"brech",
		"erschreck",
		"erschrock",
		"helf",
		"nehm",
		"schelt",
		"sprech",
		"stech",
		"sterb",
		"treff",
		"verderb",
		"verdorb",
		"werb",
		"werf",
		"beweg",
		"bewog",
		"dresch",
		"fecht",
		"flecht",
		"heb",
		"melk",
		"pfleg",
		"quell",
		"scher",
		"schmelz",
		"schwell",
		"web",
		"steh",
		"bieg",
		"biet",
		"flieg",
		"flieh",
		"fließ",
		"frier",
		"gieß",
		"kriech",
		"riech",
		"schieb",
		"schieß",
		"schließ",
		"sied",
		"sprieß",
		"stieb",
		"trief",
		"verdrieß",
		"verdross",
		"verlier",
		"verlor",
		"wieg",
		"zieh",
		"lieg",
		"beginn",
		"begonn",
		"schwimm",
		"rinn",
		"sinn",
		"spinn",
		"glimm",
		"klimm",
		"bind",
		"ding",
		"dring",
		"find",
		"kling",
		"ring",
		"schling",
		"schwind",
		"schwing",
		"sing",
		"sink",
		"spring",
		"stink",
		"trink",
		"wind",
		"wring",
		"zwing",
		"bitt",
		"sitz",
		"schind",
		"bleib",
		"leih",
		"meid",
		"preis",
		"reib",
		"scheid",
		"schein",
		"schi",
		"schreib",
		"schrei",
		"schweig",
		"spei",
		"steig",
		"treib",
		"weis",
		"verzeih",
		"verzieh",
		"beiß",
		"bleich",
		"gleich",
		"gleit",
		"greif",
		"kneif",
		"leid",
		"pfeif",
		"reiß",
		"reit",
		"scheiß",
		"schleich",
		"schleif",
		"schmeiß",
		"schneid",
		"schreit",
		"streich",
		"streit",
		"weich",
		"heiß",
		"sauf",
		"saug",
		"schnaub",
		"hau",
		"lauf",
		"komm",
		"stoß",
		"ruf",
		"häng",
		"erlösch",
		"erlosch",
		"schwör",
		"lüg",
		"trüg",
		"mahl",
		"salz",
		"spalt",
		"wiss",
		"brenn",
		"bring",
		"denk",
		"kenn",
		"nenn",
		"renn",
		"send",
		"wend",
		"hab",
		"werd",
		"könn",
		"mög",
		"dürf",
		"müss",
		"soll",
		"woll",
	
	];
	
	if(this.installedVerb.length > 0) {
		var s = this.installedVerb;
		this.installedVerb = "";
		return s;
	} else {
		return words[this.rnd() % words.length];
	}
}
	
Tarot.prototype.substantiv = function() {
	var words = [
		["","er","Ärger","en","en"],
		["","er","Abend","en","en"],
		["","er","Abflug","en","en"],
		["","er","Absender","en","en"],
		["","er","Alkohol","en","en"],
		["","er","Anfang","en","en"],
		["","er","Angestellte","en","en"],
		["","er","Anruf","en","en"],
		["","er","Anrufbeantworter","en","en"],
		["","er","Anschluss","en","en"],
		["","er","Anzug","en","en"],
		["","er","Apfel","en","en"],
		["","er","Appetit","en","en"],
		["","er","April","en","en"],
		["","er","Arbeitsplatz","en","en"],
		["","er","Arm","en","en"],
		["","er","Arzt","en","en"],
		["","er","Aufenthalt","en","en"],
		["","er","Aufzug","en","en"],
		["","er","August","en","en"],
		["","er","Ausflug","en","en"],
		["","er","Ausgang","en","en"],
		["","er","Ausländer","en","en"],
		["","er","Ausweis","en","en"],
		["","er","Automat","en","en"],
		["","er","Bahnhof","en","en"],
		["","er","Bahnsteig","en","en"],
		["","er","Balkon","en","en"],
		["","er","Baum","en","en"],
		["","er","Beamte","en","en"],
		["","er","Bekannte","en","en"],
		["","er","Berg","en","en"],
		["","er","Beruf","en","en"],
		["","er","Besuch","en","en"],
		["","er","Betrag","en","en"],
		["","er","Bildschirm","en","en"],
		["","er","Bleistift","en","en"],
		["","er","Blick","en","en"],
		["","er","Bogen","en","en"],
		["","er","Brief","en","en"],
		["","er","Briefkasten","en","en"],
		["","er","Briefumschlag","en","en"],
		["","er","Bruder","en","en"],
		["","er","Buchstabe","en","en"],
		["","er","Bus","en","en"],
		["","er","Chef","en","en"],
		["","er","Computer","en","en"],
		["","er","Dank","en","en"],
		["","er","Deutsche","en","en"],
		["","er","Dezember","en","en"],
		["","er","Dienstag","en","en"],
		["","er","Doktor","en","en"],
		["","er","Dom","en","en"],
		["","er","Donnerstag","en","en"],
		["","er","Drucker","en","en"],
		["","er","Durst","en","en"],
		["","er","Ehemann","en","en"],
		["","er","Eingang","en","en"],
		["","er","Eintritt","en","en"],
		["","er","Einwohner","en","en"],
		["","er","Empfänger","en","en"],
		["","er","Empfang","en","en"],
		["","er","Enkel","en","en"],
		["","er","Erwachsene","en","en"],
		["","er","Export","en","en"],
		["","er","Führerschein","en","en"],
		["","er","Fahrer","en","en"],
		["","er","Fahrplan","en","en"],
		["","er","Familienname","en","en"],
		["","er","Familienstand","en","en"],
		["","er","Februar","en","en"],
		["","er","Fehler","en","en"],
		["","er","Film","en","en"],
		["","er","Fisch","en","en"],
		["","er","Flughafen","en","en"],
		["","er","Flur","en","en"],
		["","er","Fluss","en","en"],
		["","er","Fotoapparat","en","en"],
		["","er","Frühling","en","en"],
		["","er","Freitag","en","en"],
		["","er","Freund","en","en"],
		["","er","Friseur","en","en"],
		["","er","Fuß","en","en"],
		["","er","Fußball","en","en"],
		["","er","Garten","en","en"],
		["","er","Gast","en","en"],
		["","er","Geburtsort","en","en"],
		["","er","Geburtstag","en","en"],
		["","er","Glückwunsch","en","en"],
		["","er","Großvater","en","en"],
		["","er","Gruß","en","en"],
		["","er","Hals","en","en"],
		["","er","Haushalt","en","en"],
		["","er","Hausmann","en","en"],
		["","er","Herbst","en","en"],
		["","er","Herd","en","en"],
		["","er","Herr","en","en"],
		["","er","Hund","en","en"],
		["","er","Hunger","en","en"],
		["","er","Import","en","en"],
		["","er","Inhalt","en","en"],
		["","er","Januar","en","en"],
		["","er","Job","en","en"],
		["","er","Jugendliche","en","en"],
		["","er","Juli","en","en"],
		["","er","Junge","en","en"],
		["","er","Juni","en","en"],
		["","er","Käse","en","en"],
		["","er","Körper","en","en"],
		["","er","Kühlschrank","en","en"],
		["","er","Kaffee","en","en"],
		["","er","Kalender","en","en"],
		["","er","Kassettenrecorder","en","en"],
		["","er","Keller","en","en"],
		["","er","Kellner","en","en"],
		["","er","Kindergarten","en","en"],
		["","er","Kinderwagen","en","en"],
		["","er","Kiosk","en","en"],
		["","er","Koffer","en","en"],
		["","er","Kollege","en","en"],
		["","er","Kontakt","en","en"],
		["","er","Kopf","en","en"],
		["","er","Kredit","en","en"],
		["","er","Kreis","en","en"],
		["","er","Kuchen","en","en"],
		["","er","Kugelschreiber","en","en"],
		["","er","Kunde","en","en"],
		["","er","Kurs","en","en"],
		["","er","Löffel","en","en"],
		["","er","Laden","en","en"],
		["","er","Lehrer","en","en"],
		["","er","LKW","en","en"],
		["","er","Loch","en","en"],
		["","er","Lohn","en","en"],
		["","er","März","en","en"],
		["","er","Müll","en","en"],
		["","er","Magen","en","en"],
		["","er","Mai","en","en"],
		["","er","Mann","en","en"],
		["","er","Mantel","en","en"],
		["","er","Markt","en","en"],
		["","er","Mechaniker","en","en"],
		["","er","Mensch","en","en"],
		["","er","Messer","en","en"],
		["","er","Mittag","en","en"],
		["","er","Mittwoch","en","en"],
		["","er","Moment","en","en"],
		["","er","Monat","en","en"],
		["","er","Montag","en","en"],
		["","er","Morgen","en","en"],
		["","er","Motor","en","en"],
		["","er","Mund","en","en"],
		["","er","Nachbar","en","en"],
		["","er","Nachmittag","en","en"],
		["","er","Name","en","en"],
		["","er","Nebel","en","en"],
		["","er","Norden","en","en"],
		["","er","Notarzt","en","en"],
		["","er","Notfall","en","en"],
		["","er","November","en","en"],
		["","er","Ober","en","en"],
		["","er","Oktober","en","en"],
		["","er","Opa","en","en"],
		["","er","Ort","en","en"],
		["","er","Osten","en","en"],
		["","er","Park","en","en"],
		["","er","Partner","en","en"],
		["","er","Pass","en","en"],
		["","er","Pkw","en","en"],
		["","er","Plan","en","en"],
		["","er","Platz","en","en"],
		["","er","Preis","en","en"],
		["","er","Prospekt","en","en"],
		["","er","Pullover","en","en"],
		["","er","Rücken","en","en"],
		["","er","Rabatt","en","en"],
		["","er","Raucher","en","en"],
		["","er","Raum","en","en"],
		["","er","Regen","en","en"],
		["","er","Reifen","en","en"],
		["","er","Reis","en","en"],
		["","er","Reiseführer","en","en"],
		["","er","Rock","en","en"],
		["","er","Rosé","en","en"],
		["","er","Rundgang","en","en"],
		["","er","Süden","en","en"],
		["","er","Saft","en","en"],
		["","er","Salat","en","en"],
		["","er","Satz","en","en"],
		["","er","Schüler","en","en"],
		["","er","Schalter","en","en"],
		["","er","Schild","en","en"],
		["","er","Schinken","en","en"],
		["","er","Schirm","en","en"],
		["","er","Schlüssel","en","en"],
		["","er","Schluss","en","en"],
		["","er","Schnee","en","en"],
		["","er","Schnupfen","en","en"],
		["","er","Schrank","en","en"],
		["","er","Schuh","en","en"],
		["","er","See","en","en"],
		["","er","September","en","en"],
		["","er","Service","en","en"],
		["","er","Sessel","en","en"],
		["","er","Sohn","en","en"],
		["","er","Sommer","en","en"],
		["","er","Sonntag","en","en"],
		["","er","Spaß","en","en"],
		["","er","Spaziergang","en","en"],
		["","er","Spielplatz","en","en"],
		["","er","Stempel","en","en"],
		["","er","Stock","en","en"],
		["","er","Stoff","en","en"],
		["","er","Strand","en","en"],
		["","er","Strom","en","en"],
		["","er","Student","en","en"],
		["","er","Stuhl","en","en"],
		["","er","Supermarkt","en","en"],
		["","er","Tag","en","en"],
		["","er","Teil","en","en"],
		["","er","Teller","en","en"],
		["","er","Teppich","en","en"],
		["","er","Termin","en","en"],
		["","er","Test","en","en"],
		["","er","Text","en","en"],
		["","er","Tipp","en","en"],
		["","er","Tisch","en","en"],
		["","er","Topf","en","en"],
		["","er","Tourist","en","en"],
		["","er","Turm","en","en"],
		["","er","Tyr","en","en"],
		["","er","Unfall","en","en"],
		["","er","Unterricht","en","en"],
		["","er","Unterschied","en","en"],
		["","er","Urlaub","en","en"],
		["","er","Vater","en","en"],
		["","er","Verein","en","en"],
		["","er","Verkäufer","en","en"],
		["","er","Verkehr","en","en"],
		["","er","Vermieter","en","en"],
		["","er","Vertrag","en","en"],
		["","er","Vogel","en","en"],
		["","er","Vormittag","en","en"],
		["","er","Vorname","en","en"],
		["","er","Wagen","en","en"],
		["","er","Wald","en","en"],
		["","er","Weg","en","en"],
		["","er","Wein","en","en"],
		["","er","Westen","en","en"],
		["","er","Wetter","en","en"],
		["","er","Wind","en","en"],
		["","er","Winter","en","en"],
		["","er","Wochentag","en","en"],
		["","er","Wunsch","en","en"],
		["","er","Zahn","en","en"],
		["","er","Zettel","en","en"],
		["","er","Zoll","en","en"],
		["","er","Zucker","en","en"],
		["","er","Zug","en","en"],





		["e","e","Ärztin","e","e"],
		["e","e","Übernachtung","e","e"],
		["e","e","Abfahrt","e","e"],
		["e","e","Adresse","e","e"],
		["e","e","Ampel","e","e"],
		["e","e","Angestellte","e","e"],
		["e","e","Angst","e","e"],
		["e","e","Ankunft","e","e"],
		["e","e","Anmeldung","e","e"],
		["e","e","Anrede","e","e"],
		["e","e","Ansage","e","e"],
		["e","e","Antwort","e","e"],
		["e","e","Anzeige","e","e"],
		["e","e","Apotheke","e","e"],
		["e","e","Arbeit","e","e"],
		["e","e","Aufgabe","e","e"],
		["e","e","Auge","e","e"],
		["e","e","Ausbildung","e","e"],
		["e","e","Auskunft","e","e"],
		["e","e","Ausländerin","e","e"],
		["e","e","Aussage","e","e"],
		["e","e","Ausstellung","e","e"],
		["e","e","Autobahn","e","e"],
		["e","e","Bäckerei","e","e"],
		["e","e","Bahn","e","e"],
		["e","e","Banane","e","e"],
		["e","e","Bank","e","e"],
		["e","e","Batterie","e","e"],
		["e","e","Beamte","e","e"],
		["e","e","Beamtin","e","e"],
		["e","e","Bekannte","e","e"],
		["e","e","Beratung","e","e"],
		["e","e","Berufsschule","e","e"],
		["e","e","Bewerbung","e","e"],
		["e","e","Birne","e","e"],
		["e","e","Bitte","e","e"],
		["e","e","Blume","e","e"],
		["e","e","Bluse","e","e"],
		["e","e","Bohne","e","e"],
		["e","e","Brücke","e","e"],
		["e","e","Briefmarke","e","e"],
		["e","e","Brieftasche","e","e"],
		["e","e","Brille","e","e"],
		["e","e","Butter","e","e"],
		["e","e","CD","e","e"],
		["e","e","CD","e","e"],
		["e","e","Creme","e","e"],
		["e","e","Dame","e","e"],
		["e","e","Dauer","e","e"],
		["e","e","Deutsche","e","e"],
		["e","e","Disco","e","e"],
		["e","e","Durchsage","e","e"],
		["e","e","Dusche","e","e"],
		["e","e","Ecke","e","e"],
		["e","e","Ehefrau","e","e"],
		["e","e","Einführung","e","e"],
		["e","e","Einladung","e","e"],
		["e","e","Eltern","e","e"],
		["e","e","Entschuldigung","e","e"],
		["e","e","Erfahrung","e","e"],
		["e","e","Erlaubnis","e","e"],
		["e","e","Ermäßigung","e","e"],
		["e","e","Erwachsene","e","e"],
		["e","e","Fähre","e","e"],
		["e","e","Führung","e","e"],
		["e","e","Fabrik","e","e"],
		["e","e","Fahrkarte","e","e"],
		["e","e","Familie","e","e"],
		["e","e","Farbe","e","e"],
		["e","e","Ferien","e","e"],
		["e","e","Feuerwehr","e","e"],
		["e","e","Firma","e","e"],
		["e","e","Flasche","e","e"],
		["e","e","Flur","e","e"],
		["e","e","Frage","e","e"],
		["e","e","Frau","e","e"],
		["e","e","Freizeit","e","e"],
		["e","e","Freundin","e","e"],
		["e","e","Frist","e","e"],
		["e","e","Gabel","e","e"],
		["e","e","Garage","e","e"],
		["e","e","Gebühr","e","e"],
		["e","e","Geldbörse","e","e"],
		["e","e","Gesamtschule","e","e"],
		["e","e","Geschwister","e","e"],
		["e","e","Gesundheit","e","e"],
		["e","e","Größe","e","e"],
		["e","e","Grippe","e","e"],
		["e","e","Großeltern","e","e"],
		["e","e","Großmutter","e","e"],
		["e","e","Grundschule","e","e"],
		["e","e","Gruppe","e","e"],
		["e","e","Halbpension","e","e"],
		["e","e","Halle","e","e"],
		["e","e","Haltestelle","e","e"],
		["e","e","Hand","e","e"],
		["e","e","Hausaufgabe","e","e"],
		["e","e","Hausfrau","e","e"],
		["e","e","Heimat","e","e"],
		["e","e","Heizung","e","e"],
		["e","e","Hilfe","e","e"],
		["e","e","Hose","e","e"],
		["e","e","Idee","e","e"],
		["e","e","Industrie","e","e"],
		["e","e","Information","e","e"],
		["e","e","Jacke","e","e"],
		["e","e","Jugendherberge","e","e"],
		["e","e","Jugendliche","e","e"],
		["e","e","Küche","e","e"],
		["e","e","Kündigung","e","e"],
		["e","e","Kamera","e","e"],
		["e","e","Kanne","e","e"],
		["e","e","Karte","e","e"],
		["e","e","Kartoffel","e","e"],
		["e","e","Kasse","e","e"],
		["e","e","Kassette","e","e"],
		["e","e","Katze","e","e"],
		["e","e","Kenntnisse","e","e"],
		["e","e","Kette","e","e"],
		["e","e","Kirche","e","e"],
		["e","e","Klasse","e","e"],
		["e","e","Kleidung","e","e"],
		["e","e","Kneipe","e","e"],
		["e","e","Kollegin","e","e"],
		["e","e","Kontrolle","e","e"],
		["e","e","Kosmetik","e","e"],
		["e","e","Krankenkasse","e","e"],
		["e","e","Krankheit","e","e"],
		["e","e","Kreditkarte","e","e"],
		["e","e","Kreuzung","e","e"],
		["e","e","Kunde","e","e"],
		["e","e","Kundin","e","e"],
		["e","e","Lösung","e","e"],
		["e","e","Lampe","e","e"],
		["e","e","Landschaft","e","e"],
		["e","e","Lebensmittel","e","e"],
		["e","e","Lehre","e","e"],
		["e","e","Lehrerin","e","e"],
		["e","e","Leute","e","e"],
		["e","e","Lohn","e","e"],
		["e","e","Luft","e","e"],
		["e","e","Lust","e","e"],
		["e","e","Möbel","e","e"],
		["e","e","Mülltonne","e","e"],
		["e","e","Maschine","e","e"],
		["e","e","Mehrwertsteuer","e","e"],
		["e","e","Meinung","e","e"],
		["e","e","Menge","e","e"],
		["e","e","Miete","e","e"],
		["e","e","Milch","e","e"],
		["e","e","Minute","e","e"],
		["e","e","Mitte","e","e"],
		["e","e","Mitteilung","e","e"],
		["e","e","Mittelschule","e","e"],
		["e","e","Mode","e","e"],
		["e","e","Musik","e","e"],
		["e","e","Mutter","e","e"],
		["e","e","Nähe","e","e"],
		["e","e","Nachbarin","e","e"],
		["e","e","Nacht","e","e"],
		["e","e","Natur","e","e"],
		["e","e","Note","e","e"],
		["e","e","Notiz","e","e"],
		["e","e","Nudel","e","e"],
		["e","e","Nummer","e","e"],
		["e","e","Oma","e","e"],
		["e","e","Operation","e","e"],
		["e","e","Orange","e","e"],
		["e","e","Ordnung","e","e"],
		["e","e","Panne","e","e"],
		["e","e","Partei","e","e"],
		["e","e","Partnerin","e","e"],
		["e","e","Party","e","e"],
		["e","e","Pause","e","e"],
		["e","e","Pension","e","e"],
		["e","e","Plastik","e","e"],
		["e","e","Polizei","e","e"],
		["e","e","Portion","e","e"],
		["e","e","Post","e","e"],
		["e","e","Postleitzahl","e","e"],
		["e","e","Prüfung","e","e"],
		["e","e","Praxis","e","e"],
		["e","e","Qualität","e","e"],
		["e","e","Quittung","e","e"],
		["e","e","Raucherin","e","e"],
		["e","e","Realschule","e","e"],
		["e","e","Rechnung","e","e"],
		["e","e","Reinigung","e","e"],
		["e","e","Reise","e","e"],
		["e","e","Reparatur","e","e"],
		["e","e","Rezeption","e","e"],
		["e","e","Rosé","e","e"],
		["e","e","Rose","e","e"],
		["e","e","Sache","e","e"],
		["e","e","Schülerin","e","e"],
		["e","e","Scheckkarte","e","e"],
		["e","e","Schokolade","e","e"],
		["e","e","Schule","e","e"],
		["e","e","Schwester","e","e"],
		["e","e","See","e","e"],
		["e","e","Sehenswürdigkeit","e","e"],
		["e","e","Seife","e","e"],
		["e","e","Sekretärin","e","e"],
		["e","e","Sekunde","e","e"],
		["e","e","Sendung","e","e"],
		["e","e","Sonne","e","e"],
		["e","e","Sorge","e","e"],
		["e","e","Spülmaschine","e","e"],
		["e","e","Speisekarte","e","e"],
		["e","e","Sprache","e","e"],
		["e","e","Sprachschule","e","e"],
		["e","e","Sprechstunde","e","e"],
		["e","e","Stadt","e","e"],
		["e","e","Steuer","e","e"],
		["e","e","Straße","e","e"],
		["e","e","Straßenbahn","e","e"],
		["e","e","Studentin","e","e"],
		["e","e","Stunde","e","e"],
		["e","e","Suppe","e","e"],
		["e","e","Tür","e","e"],
		["e","e","Tüte","e","e"],
		["e","e","Tankstelle","e","e"],
		["e","e","Tasche","e","e"],
		["e","e","Tasse","e","e"],
		["e","e","Tochter","e","e"],
		["e","e","Toilette","e","e"],
		["e","e","Tomate","e","e"],
		["e","e","Treppe","e","e"],
		["e","e","Tür","e","e"],
		["e","e","Uhr","e","e"],
		["e","e","Universität","e","e"],
		["e","e","Unterhaltung","e","e"],
		["e","e","Unterkunft","e","e"],
		["e","e","Unterschrift","e","e"],
		["e","e","Untersuchung","e","e"],
		["e","e","Verbindung","e","e"],
		["e","e","Verkäuferin","e","e"],
		["e","e","Versicherung","e","e"],
		["e","e","Verspätung","e","e"],
		["e","e","Volksschule","e","e"],
		["e","e","Vorsicht","e","e"],
		["e","e","Vorwahl","e","e"],
		["e","e","Wäsche","e","e"],
		["e","e","Welt","e","e"],
		["e","e","Werkstatt","e","e"],
		["e","e","Wirtschaft","e","e"],
		["e","e","Woche","e","e"],
		["e","e","Wohnung","e","e"],
		["e","e","Wolke","e","e"],
		["e","e","Wurst","e","e"],
		["e","e","Zahl","e","e"],
		["e","e","Zeit","e","e"],
		["e","e","Zeitschrift","e","e"],
		["e","e","Zeitung","e","e"],
		["e","e","Zigarette","e","e"],
		["e","e","Zitrone","e","e"],



		["","es","Alter","","es"],
		["","es","Angebot","","es"],
		["","es","Appartement","","es"],
		["","es","Auge","","es"],
		["","es","Ausland","","es"],
		["","es","Auto","","es"],
		["","es","Büro","","es"],
		["","es","Baby","","es"],
		["","es","Bad","","es"],
		["","es","Bein","","es"],
		["","es","Beispiel","","es"],
		["","es","Benzin","","es"],
		["","es","Bett","","es"],
		["","es","Bier","","es"],
		["","es","Bild","","es"],
		["","es","Blatt","","es"],
		["","es","Blut","","es"],
		["","es","Brötchen","","es"],
		["","es","Brot","","es"],
		["","es","Buch","","es"],
		["","es","Cafe","","es"],
		["","es","Dach","","es"],
		["","es","Datum","","es"],
		["","es","Ding","","es"],
		["","es","Doppelzimmer","","es"],
		["","es","Dorf","","es"],
		["","es","Ei","","es"],
		["","es","Einzelzimmer","","es"],
		["","es","Eis","","es"],
		["","es","Ende","","es"],
		["","es","Erdgeschoss","","es"],
		["","es","Ergebnis","","es"],
		["","es","Essen","","es"],
		["","es","Fahrrad","","es"],
		["","es","Fax","","es"],
		["","es","Fenster","","es"],
		["","es","Fernsehgerät","","es"],
		["","es","Fest","","es"],
		["","es","Feuer","","es"],
		["","es","Feuerzeug","","es"],
		["","es","Fieber","","es"],
		["","es","Fleisch","","es"],
		["","es","Flugzeug","","es"],
		["","es","Formular","","es"],
		["","es","Foto","","es"],
		["","es","Frühjahr","","es"],
		["","es","Frühstück","","es"],
		["","es","Fundbüro","","es"],
		["","es","Gas","","es"],
		["","es","Geburtsjahr","","es"],
		["","es","Gegenteil","","es"],
		["","es","Geld","","es"],
		["","es","Gemüse","","es"],
		["","es","Gepäck","","es"],
		["","es","Gericht","","es"],
		["","es","Geschäft","","es"],
		["","es","Geschenk","","es"],
		["","es","Geschirr","","es"],
		["","es","Gesicht","","es"],
		["","es","Gespräch","","es"],
		["","es","Getränk","","es"],
		["","es","Gewicht","","es"],
		["","es","Gewitter","","es"],
		["","es","Glück","","es"],
		["","es","Glas","","es"],
		["","es","Gleis","","es"],
		["","es","Guthaben","","es"],
		["","es","Gymnasium","","es"],
		["","es","Hähnchen","","es"],
		["","es","Haar","","es"],
		["","es","Handtuch","","es"],
		["","es","Handy","","es"],
		["","es","Haus","","es"],
		["","es","Hemd","","es"],
		["","es","Herz","","es"],
		["","es","Hobby","","es"],
		["","es","Holz","","es"],
		["","es","Internet","","es"],
		["","es","Jahr","","es"],
		["","es","Junge","","es"],
		["","es","Kennzeichen","","es"],
		["","es","Kfz","","es"],
		["","es","Kind","","es"],
		["","es","Kino","","es"],
		["","es","Kleid","","es"],
		["","es","Konsulat","","es"],
		["","es","Konto","","es"],
		["","es","Konzert","","es"],
		["","es","Laden","","es"],
		["","es","Lager","","es"],
		["","es","Land","","es"],
		["","es","Leben","","es"],
		["","es","Lebensmittel","","es"],
		["","es","Leid","","es"],
		["","es","Licht","","es"],
		["","es","Lied","","es"],
		["","es","Loch","","es"],
		["","es","Lokal","","es"],
		["","es","Mädchen","","es"],
		["","es","Möbel","","es"],
		["","es","Mal","","es"],
		["","es","Material","","es"],
		["","es","Medikament","","es"],
		["","es","Meer","","es"],
		["","es","Messer","","es"],
		["","es","Metall","","es"],
		["","es","Mittel","","es"],
		["","es","Moment","","es"],
		["","es","Morgen","","es"],
		["","es","Museum","","es"],
		["","es","Obst","","es"],
		["","es","Öl","","es"],
		["","es","Päckchen","","es"],
		["","es","Paket","","es"],
		["","es","Papier","","es"],
		["","es","Parfüm","","es"],
		["","es","Plastik","","es"],
		["","es","Praktikum","","es"],
		["","es","Problem","","es"],
		["","es","Programm","","es"],
		["","es","Rücken","","es"],
		["","es","Radio","","es"],
		["","es","Rathaus","","es"],
		["","es","Reifen","","es"],
		["","es","Reis","","es"],
		["","es","Reisebüro","","es"],
		["","es","Restaurant","","es"],
		["","es","Rezept","","es"],
		["","es","Rind","","es"],
		["","es","Rosé","","es"],
		["","es","Salz","","es"],
		["","es","Schiff","","es"],
		["","es","Schild","","es"],
		["","es","Schloss","","es"],
		["","es","Schwein","","es"],
		["","es","Schwimmbad","","es"],
		["","es","Service","","es"],
		["","es","Sofa","","es"],
		["","es","Sonderangebot","","es"],
		["","es","Stück","","es"],
		["","es","Standesamt","","es"],
		["","es","Steuer","","es"],
		["","es","Streichholz","","es"],
		["","es","Studium","","es"],
		["","es","Tag","","es"],
		["","es","Taxi","","es"],
		["","es","Teil","","es"],
		["","es","Telefon","","es"],
		["","es","Telefonbuch","","es"],
		["","es","Theater","","es"],
		["","es","Thema","","es"],
		["","es","Ticket","","es"],
		["","es","Tier","","es"],
		["","es","Trinkgeld","","es"],
		["","es","Video","","es"],
		["","es","Wasser","","es"],
		["","es","Werkzeug","","es"],
		["","es","Wetter","","es"],
		["","es","Wiederhören","","es"],
		["","es","Wiedersehen","","es"],
		["","es","Wochenende","","es"],
		["","es","Wort","","es"],
		["","es","Zentrum","","es"],
		["","es","Zeugnis","","es"],
		["","es","Zimmer","","es"],





	
		["","er","Baum","en","en"], 
		["e","e","Sympathie","e","e"], 
		["e","e","Zahl","e","e"], 
		["","er","Umlauf","en","en"], 
		["","er","Stopp","","es"], 
		["e","e","Liebe","e","e"],
		["e","e","Decke","e","e"],
		["e","e","Wand","e","e"],
		["e","e","Woche","e","e"],
		["e","e","Blume","e","e"],
		["","er","Hass","en","en"], 
		["","er","Himmel","en","en"], 
		["","er","Pfleger","en","en"], 
		["e","e","Hölle","e","e"],
		["e","e","Pflegerin","e","e"],
		["e","e","Tastatur","e","e"],
		["e","e","Betreuerin","e","e"],
		["e","e","Arbeitgeberin","e","e"],
		["e","e","Arbeitnehmerin","e","e"],
		["e","e","Betreuerin","e","e"],
		["","er","Monat","en","en"], 
		["","er","Tag","en","en"], 
		["","er","Herr Gott","en","en"], 
		["","er","Klient","en","en"], 
		["","er","Arbeitgeber","en","en"], 
		["","er","Arbeitnehmer","en","en"], 
		["","er","Betreuer","en","en"], 
		["","er","Mund","en","en"], 
		["","er","Gott","en","en"], 
		["","er","Garten","en","en"], 
		["","es","Lamm","","es"],
		["","es","Blümchen","","es"],
		["","es","Land","","es"],
		["","es","Haus","","es"],
		["","es","Jahr","","es"],
		["","es","Dach","","es"],
		["","er","Teufel","en","en"], 
		["","er","Computer","en","en"], 
		["","es","Geld","","es"], 
		["","er","Boden","en","en"], 
		["","er","Zufall","en","en"], 
		["","er","Kopf","en","en"], 
		["","er","Bauch","en","en"], 
		["","er","Hunger","en","en"], 
		["","er","Finger","en","en"], 
		["","er","Fuß","en","en"], 
		["","er","Tee","en","en"], 
		["","er","Laden","en","en"], 
		["e","e","Informatik","e","e"],
		["e","e","Ironie","e","e"],
		["e","e","Praxis","e","e"],
		["e","e","Theorie","e","e"],
		["e","e","Biologie","e","e"],
		["e","e","Wissenschaft","e","e"],
		["e","e","Psychologie","e","e"],
		["e","e","Partei","e","e"],
		["e","e","Justiz","e","e"],
		["e","e","Regierung","e","e"],
		["e","e","Zeit","e","e"],
		["e","e","Polizei","e","e"],
		["e","e","Illusion","e","e"],
		["e","e","Arbeit","e","e"],
		["e","e","Buße","e","e"]
	];
	
	var tp = words[this.rnd() % words.length];
	var tp2 = words[this.rnd() % words.length];
	
	switch(this.rnd() % 4) {
		case 0:
			return tp;
		case 1:
			return [tp2[0],tp2[1],tp[2] + "-" + tp2[2],tp2[3],tp2[4]];
		case 2:
			return ["e","e",this.upFirst(this.verb() + "ung"),"e","e"];
		case 3:
			return [tp[0],tp[1],this.verb() + "-" + tp[2],tp[3],tp[4]];
	}
}
	
Tarot.prototype.upFirst = function(s) {
	return s.substr(0,1).toUpperCase() + s.substr(1);
}
	
Tarot.prototype.adjektiv = function() {
	var words = [
		"abstinent",
		"achtsam",
		"afrikanisch",
		"akkurat",
		"alkoholisch",
		"alphabetisch",
		"ängstlich",
		"ärgerlich",
		"baff",
		"bairisch",
		"bang",
		"bankrott",
		"bedrohlich",
		"begünstigt",
		"behaglich",
		"beharrlich",
		"blind",
		"brillant",
		"charmant",
		"chemisch",
		"chorisch",
		"christlich",
		"chronisch",
		"chronologisch",
		"dämlich",
		"dämmrig",
		"dankbar",
		"darstellbar",
		"dazugehörig",
		"deckend",
		"demokratisch",
		"depressiv",
		"derb",
		"dialogisch",
		"diebisch",
		"dumm",
		"dringend",
		"eckig",
		"edel",
		"effizient",
		"egoistisch",
		"ehebrecherisch",
		"ehrerbietig",
		"ehrfürchtig",
		"ehrgeizig",
		"ehrlos",
		"eigenständig",
		"einladend",
		"elektrisch",
		"fabelhaft",
		"fachkundig",
		"fad", 
		"fadenscheinig",
		"fahrlässig",
		"faktisch",
		"fantasielos",
		"fantastisch",
		"fein",
		"fest",
		"fettig",
		"fit",
		"flach",
		"flauschig",
		"flott",
		"gängig",
		"garstig",
		"gastfreundlich",
		"gebogen",
		"gedrückt",
		"geeignet",
		"gefährlich",
		"gefangen",
		"geisterhaft",
		"gelb",
		"gereizt",
		"glatt",
		"gleichberechtigt",
		"glücklich",
		"grafisch",
		"groß",
		"haarig",
		"halb",
		"halsbrecherisch",
		"hämisch",
		"handlungsfähig",
		"heiß",
		"hell",
		"herzoglich",
		"hinfällig",
		"hoch",
		"hoffnungsvoll",
		"hündisch",
		"ideal",
		"identisch",
		"idyllisch",
		"illegal",
		"illusorisch",
		"imaginär",
		"imponierend",
		"individuell",
		"inhaltlich",
		"inklusiv",
		"integriert",
		"international",
		"jährlich",
		"jetzig",
		"jodhaltig",
		"jordanisch",
		"jüdisch",
		"jugendlich",
		"jung",
		"jungfräulich",
		"kahl",
		"kalorisch",
		"kämpferisch",
		"katholisch",
		"käuflich",
		"keusch",
		"kirchlich",
		"klangvoll",
		"knackig",
		"kokett",
		"kontrovers",
		"korrekt",
		"krank",
		"krumm",
		"künstlich",
		"labberig", 
		"labil",
		"lahm",
		"ländlich",
		"laut",
		"lebensgroß",
		"legitim",
		"leicht",
		"lieb",
		"lockig",
		"lokal",
		"löslich",
		"luftig",
		"mächtig",
		"männlich",
		"maßvoll",
		"materiell",
		"mehrsprachig",
		"meisterlich",
		"mental",
		"mickerig", 
		"mitleidig",
		"monatlich",
		"motorisch",
		"musikalisch",
		"mutig",
		"mütterlich",
		"nächtlich",
		"nah",
		"närrisch",
		"nass",
		"negativ",
		"neidisch",
		"neu",
		"niedrig",
		"niveaulos",
		"nördlich",
		"normal",
		"notdürftig",
		"oberschlau",
		"obsolet",
		"obszön",
		"ockerfarben", 
		"öde", 
		"offen",
		"öffentlich",
		"ökologisch",
		"ölig",
		"olympiareif",
		"operativ",
		"oral",
		"örtlich",
		"österlich",
		"pädagogisch",
		"paradiesisch",
		"parkartig",
		"parlamentarisch",
		"passiv",
		"peinlich",
		"pensioniert",
		"persönlich",
		"perspektivlos",
		"pflichtbewusst",
		"phantastisch",
		"physisch",
		"politisch",
		"poetisch",
		"praktisch",
		"provokant",
		"qualitativ",
		"qualvoll",
		"quecksilberhaltig",
		"quengelig",
		"quergestreift",
		"quicklebendig",
		"quietschfidel",
		"rabenschwarz",
		"radikal",
		"raffiniert",
		"rankenartig",
		"rasch",
		"ratlos",
		"rauchfrei",
		"recyclebar",
		"reformbedürftig",
		"regnerisch",
		"reich",
		"rein",
		"relativ",
		"respektvoll", 
		"rhythmisch",
		"riesig",
		"roh",
		"rostig",
		"rückläufig",
		"sachkundig",
		"sachlich",
		"saisonal",
		"salzig",
		"sauer",
		"scharf",
		"schattig",
		"schleimig",
		"schreckenerregend",
		"schusselig",
		"seidenweich",
		"selbstständig",
		"sesshaft",
		"sicher",
		"soft",
		"sperrig",
		"spitz",
		"steil",
		"stramm",
		"tailliert",
		"taktvoll",
		"technisch",
		"temperamentvoll",
		"theoretisch",
		"topografisch",
		"tot",
		"trächtig",
		"traditionell",
		"treu",
		"trocken",
		"trotzig",
		"tüchtig",
		"übel",
		"übertrieben",
		"überparteilich",
		"ultimativ",
		"ultrakurz",
		"umkehrbar",
		"umständlich",
		"unbehaglich",
		"unerlässlich",
		"uralt",
		"variabel",
		"väterlich",
		"verabscheuungswürdig",
		"verantwortungslos",
		"verblüfft",
		"verdaulich",
		"verklemmt",
		"versichert",
		"viertürig",
		"viral",
		"voll",
		"wach",
		"weise",
		"wahlberechtigt",
		"warm",
		"wässrig",
		"weich",
		"weihnachtlich",
		"weit",
		"weise",
		"weiß",
		"weitreichend",
		"wertvoll",
		"widerlegbar",
		"wirtschaftlich",
		"wohnlich",
		"x,-beliebig",
		"x-fach",
		"y-förmig",
		"zackig",
		"zahlenmäßig",
		"zappelig", 
		"zart",
		"zehnjährig",
		"zeitlich",
		"zentral",
		"zickig",
		"zinslos",
		"zivil",
		"zornig",
		"zuckerfrei",
		"zuvorkommend",
		"zweckgebunden",
		"zweifach",
		"zynisch",		
	];
	
	return words[this.rnd() % words.length];
}
	
Tarot.prototype.attribut = function(adj) {
	if(adj.length == 0)
		switch(this.rnd() % 2) {
			case 0:
				return this.adjektiv();
			case 1:
				return "ge" + this.verb() + "et";
		}
	else
		return adj;
}
	
Tarot.prototype.frageWort = function() {
	var words = [
		"Für wen", "Wie", "Für was", "Warum", "Wann", "Wo"
	];
	
	return words[this.rnd() % words.length];
}
	
Tarot.prototype.gebeugtesVerb = function(person, plural, tempus) {
	switch(person) {
		case 1:
			if(plural) {
				switch(tempus) {
					case -1:
						return [this.verb() + "eten",""];
					case 0:
						return [this.verb() + "en",""];
					case 1:
						return ["werden",this.verb() + "en"];
				};
			} else {
				switch(tempus) {
					case -1:
						return [this.verb() + "ete",""];
					case 0:
						return [this.verb() + "e",""];
					case 1:
						return ["werde",this.verb() + "en"];
				};
			}
		case 2:
			if(plural) {
				switch(tempus) {
					case -1:
						return ["habt", "ge" + this.verb() + "et"];
					case 0:
						return [this.verb() + "et",""];
					case 1:
						return ["werdet",this.verb() + "en"];
				};
			} else {
				switch(tempus) {
					case -1:
						return ["hast","ge" + this.verb() + "et"];
					case 0:
						return [this.verb() + "st",""];
					case 1:
						return ["wirst",this.verb() + "en"];
				};
			}
		case 3:
			if(plural) {
				switch(tempus) {
					case -1:
						return ["haben","ge" + this.verb() + "et"];
					case 0:
						return [this.verb() + "en",""];
					case 1:
						return ["werden",this.verb() + "en"];
				};
			} else {
				switch(tempus) {
					case -1:
						return ["hat","ge" + this.verb() + "et"];
					case 0:
						return [this.verb() + "et",""];
					case 1:
						return ["wird",this.verb() + "en"];
				};
			}
	};
}
	
Tarot.prototype.pronomen = function(person,plural) {
	switch(person) {
		case 1:
			if(plural) {
				return "wir";
			} else {
				return "ich";
			}
		case 2:
			if(plural) {
				return "ihr";
			} else {
				return "du";
			}
		case 3:
			if(plural) {
				return "sie";
			} else {
				var words = [
					"er", "sie", "es"
				];
				
				return words[this.rnd() % words.length];
			}
	};
}
	
Tarot.prototype.attributiert = function(plural,akkusativ,adj) {
	var subs = this.substantiv();
	
	var words = [
		"alle", "mehrere", "manche", "wenige", "viele"
	];
	var x = words[this.rnd() % words.length];
	
	var words2 = [
		"ein", "genau ein", "gerademal ein", "ungefähr ein"
	];
	var y = words2[this.rnd() % words2.length];
	
	if(akkusativ) {
		if(plural) {
			if(this.rnd() % 2) 
				return x + " " + this.attribut(adj) + "en " + subs[2] + "en";
			else
				return x + " " + this.attribut("") + " " + this.attribut(adj) + "en " + subs[2] + "en";
		} else {
			if(this.rnd() % 2) 
				return y + subs[3] + " " + this.attribut(adj) + subs[4] + " " + subs[2];
			else 
				return y + subs[3] + " " + this.attribut("") + " " + this.attribut(adj) + subs[4] + " " + subs[2];
		}
	} else {
		if(plural) {
			if(this.rnd() % 2) 
				return x + " " + this.attribut(adj) + "en " + subs[2] + "en";
			else
				return x + " " + this.attribut("") + " " + this.attribut(adj) + "en " + subs[2] + "en";
		} else {
			if(this.rnd() % 2) 
				return y + subs[0] + " " + this.attribut(adj) + subs[1] + " " + subs[2];
			else
				return y + subs[0] + " " + this.attribut("") + " " + this.attribut(adj) + subs[1] + " " + subs[2];
		}
	}
}


Tarot.prototype.zeitAngabe = function(y) {
	if(this.rnd() % 3 && !y) return "von " + this.zeitAngabe(true) + " bis " + this.zeitAngabe(true);
	
	switch(y ? this.rnd() % 3 : this.rnd() % 5) {
		case 0:
			let as = ["Montag","Dienstag","Mittwoch","Donnerstag","Freitag","Samstag","Sonntag"];
			let a = as[this.rnd() % as.length];
			
			if(this.rnd() % 2) a = this.adjektiv() + "en " + a;
				
			return "am " + a;
			
		case 1:
			let bs = ["Januar","Februar","März","April","Mai","Juni","Juli","August",
				"September","Oktober","November","Dezember"];
			let b = bs[this.rnd() % bs.length];
			
			if(this.rnd() % 2)
				b = b + " " + this.rnd() % 3000;
				
			if(this.rnd() % 2)
				b = this.rnd() % 31 + "ten " + b;

			if(this.rnd() % 2)
				b = this.adjektiv() + "en " + b;
				
			return "am " + b;
		
		case 2:
			let cs = ["gestern","vorgestern","heute","morgen","übermorgen"];
			return cs[this.rnd() % cs.length];
			
		default:
			return "";
	}
}

	
Tarot.prototype.frage = function() {
	var person = this.rnd() % 3 + 1;
	var tempus = this.rnd() % 3 - 1;
	var plural = this.rnd() % 2;
	var gv = this.gebeugtesVerb(person,plural,tempus);
	var pn = this.pronomen(person,plural);
			
	switch(person) {
		case 3:
			return this.frageWort() + " " + gv[0] + " "
				+ this.attributiert(plural,false,"") + " " 
				+ this.attributiert(this.rnd() % 2,true,"") + " " 
				+ gv[1] + " ?";
		case 1:
		case 2:
			return this.frageWort() + " " + gv[0] + " "
				+ pn + " " 
				+ this.attributiert(this.rnd() % 2,true,"") + " " 
				+ gv[1] + " ?";
	}
}
	
Tarot.prototype.antwort = function(adj1,adj2) {
	var person = this.rnd() % 3 + 1;
	if(adj1.length >= 0) person = 3;

	var tempus = this.rnd() % 3 - 1;
	var plural = this.rnd() % 2;
	var gv = this.gebeugtesVerb(person,plural,tempus);
	var pn = this.pronomen(person,plural);
	
	switch(person) {
		case 3:
			return this.zeitAngabe() + " " + this.attributiert(plural,false,adj1) + " " 
				+ gv[0] + " "
				+ this.attributiert(this.rnd() % 2,true,adj2) + " "
				+ gv[1];
		case 1:
		case 2:
			return this.zeitAngabe() + " " + pn + " " 
				+ gv[0] + " "
				+ this.attributiert(this.rnd() % 2,true,"") + " "
				+ gv[1];
	}
}
	
Tarot.prototype.joke = function() {
	var t = this.clichee();
	
	this.installedVerb = t[2];
	
	var verbs = [
		t[3] + " unter",t[3] + " über",
		"langsam " + t[3],"sehr langsam " + t[3],
		"schnell " + t[3],"sehr schnell " + t[3],
		t[3]
	];
	
	if(this.rnd() % 2) {
		var swap = t[0];
		t[0] = t[1];
		t[1] = swap;
	}
	
	return this.antwort(t[0],t[1]) + " " + verbs[this.rnd() % verbs.length] + " " + t[4];
}
	
Tarot.prototype._get = function() {
	switch(this.rnd() % 3) {
		case 0:
			return this.frage();
		case 1:
			return this.antwort("","") + ".";
		case 2:
			return this.joke() + ".";
		case 3: 
			var x = this.rnd() % 2 ? this.antwort("","") : this.joke();
			if(this.rnd() % 2)
				x += ", so dass " + this.antwort("","");
			else
				x += ", so dass " + this.joke();
			
			return x;
	}
}


Tarot.prototype.get = function() {
	return this.upFirst(this._get().trim());
}




if(typeof module !== "undefined")
module.exports = {
	md5 : md5,
	random : Random,
	tarot : Tarot,
	shuffle : shuffle,
	Chess : {
		Chess : Chess,
		ChessVariants : {
			CLassicChess : ClassicChess,
			ZauberStrikeChess : ZauberStrikeChess,
			GardenEdenChess : GardenEdenChess
		}
	},
	makeScale : makeScale,
	Pattern : Pattern
};


