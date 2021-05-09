importScripts("MyLib.js");


function myPropose(me,board,n) {
	if(n == 0)  {
		return me.measureWhite(board);
	}

	let m0 = null;

	for(let i = 0;i < 8;++i)
	for(let j = 0;j < 8;++j)
	if(board[i + 8 * j] && (board[i + 8 * j].color == board.color || me.team())) {
		let ms = me.endFilter(board,1,i,j,me.moveFilter(board,i,j,me.getMoves(board,i,j)));
		ms.forEach(mv => {
			let b = me.doMove(board,i,j,mv);
			let m = myPropose(me,b,n - 1);
			if(m0 == null) {
				m0 = m;
			}
			if(me.team()) {
				if(m > m0) {
					m0 = m;
				}

			} else {
				if(board.color == "white" && m > m0) {
					m0 = m;
				} else if(board.color == "black" && m < m0) {
					m0 = m;
				}
			}
		});
	}
	return m0;
}


onmessage = (msg) => {
	let me = null;

	switch(msg.data.type) {
		case "ClassicChess":
			me = new ClassicChess();
			break;
		case "ZauberStrikeChess":
			me = new ZauberStrikeChess();
			break;
		case "GardenEdenChess":
			me = new GardenEdenChess();
			break;
		case "HalmaChess":
			me = new HalmaChess();
			break;
		case "GalagaChess":
			me = new GalagaChess();
			break;
		case "MonopolyChess":
			me = new MonopolyChess();
			break;
	};

	me.board = msg.data.thisBoard;
	me._team = msg.data._team;
	
	let m = myPropose(me,msg.data.board,msg.data.n);
	
	msg.target.postMessage([
		msg.data.type,m
	]);
}


