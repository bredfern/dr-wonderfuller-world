----------------------------------------------------------------------------------
-- Company: 
-- Engineer: 
-- 
-- Create Date:    21:49:43 06/27/2008 
-- Design Name: 
-- Module Name:    main - Behavioral 
-- Project Name: 
-- Target Devices: 
-- Tool versions: 
-- Description: 
--
-- Dependencies: 
--
-- Revision: 
-- Revision 0.01 - File Created
-- Additional Comments: 
--
----------------------------------------------------------------------------------
library IEEE;
use IEEE.STD_LOGIC_1164.ALL;
--use IEEE.STD_LOGIC_ARITH.ALL;
--use IEEE.STD_LOGIC_UNSIGNED.ALL;
use IEEE.numeric_std.ALL;


---- Uncomment the following library declaration if instantiating
---- any Xilinx primitives in this code.
--library UNISIM;
--use UNISIM.VComponents.all;

entity main is
	port(
		clk,reset : in std_logic;
		--btn: in std_logic_vector(2 downto 0);
		--rx : in std_logic;
		--tx : out std_logic;
		led: out std_logic_vector(7 downto 0);
		--sseg : out std_logic_vector(7 downto 0);
		--an : out std_logic_vector(3 downto 0);
		sw : in std_logic_vector(2 downto 0);
		ps2d,ps2c : inout std_logic;
		
		-- to/from chip
		ad : out std_logic_vector(17 downto 0);
		we_n, oe_n : out std_logic;
		-- SRAM chip a
		dio_a : inout std_logic_vector(15 downto 0);
		ce_a_n,ub_a_n,lb_a_n : out std_logic;
		-- SRAM chip b
		dio_b : inout std_logic_vector(15 downto 0);
		ce_b_n,ub_b_n,lb_b_n : out std_logic;

		hsync,vsync : out std_logic;
		rgb : out std_logic_vector(2 downto 0)
	);
end main;

architecture Behavioral of main is
	signal pixel_tick : std_logic;
	signal wr_tick : std_logic;
	signal wr_addr : std_logic_vector(17 downto 0);
	signal wr_data : std_logic_vector(31 downto 0);
	signal wr_fifo_full : std_logic;
	signal video_on : std_logic;
	signal pixel_x,pixel_y : std_logic_vector(9 downto 0);
	signal vpu_rgb : std_logic_vector(2 downto 0);
	signal xm,ym : std_logic_vector(8 downto 0);
	signal btnm : std_logic_vector(2 downto 0);
	signal m_done_tick: std_logic;
	signal sumx,sumy : unsigned(9 downto 0);
	signal addsub_a,addsub_b,addsub_result : std_logic_vector(17 downto 0);
	signal addsub_operation : std_logic_vector(5 downto 0);
	signal addsub_operation_nd,addsub_operation_rfd,addsub_rdy : std_logic;	
	signal ftf_a :std_logic_vector(11 downto 0);
	signal ftf_operation_nd,ftf_operation_rfd,ftf_rdy : std_logic;	
	signal ftf_result : std_logic_vector(17 downto 0);
	signal mult_a,mult_b,mult_result : std_logic_vector(17 downto 0);
	signal mult_operation_nd,mult_operation_rfd,mult_rdy : std_logic;	
	signal div_a,div_b,div_result : std_logic_vector(17 downto 0);
	signal div_operation_nd,div_operation_rfd,div_rdy : std_logic;	
	signal fltf_a :std_logic_vector(17 downto 0);
	signal fltf_operation_nd,fltf_operation_rfd,fltf_rdy : std_logic;	
	signal fltf_result : std_logic_vector(11 downto 0);
	signal mouse_reset : std_logic;
	
	type state_t is (s3,s4,s5,s10,s20,s22,s24,s26,s28,s29,s30,calcxm,calcym,colorcalculation,sub_linear,initcolcal);
	type state_reset_t is (s10,s20,s30,s40,s50,s60,s70,s80,s90);
	type state_sub_linear_t is (s5,s8,s10,s20,s30,s40,s50,s60,s70,s80);
	type state_colcal_t is (s10,s20,s30,s40,s50,s60,s70,s80,s90,s100,
		s110,s120,s130,s140,s150,s160,s170,s180,s190,s200,
		s210,s220,s230,s240,s250,s260,s270);
	
	type regs_t is record
			state : state_t;
			state_reset : state_reset_t;
			state_colcal : state_colcal_t;
			state_sub_linear : state_sub_linear_t;
			address : unsigned(17 downto 0);
			lineaddress : unsigned(17 downto 0);
			n : unsigned(12 downto 0);
			x,y : unsigned(11 downto 0);
			xmax,xmin,xmath,width : std_logic_vector(17 downto 0);
			ymax,ymin,ymath,height : std_logic_vector(17 downto 0);
			r1,r2,r3,r4,r5,r6 : std_logic_vector(17 downto 0);
			return_state : state_t;
			dataword : std_logic_vector(31 downto 0);
			i : unsigned(3 downto 0);
			color : std_logic_vector(2 downto 0);
		end record;
	signal rreg,rnext : regs_t;

	type state_mouse_t is (s10,s20,s30);
	type regs2_t is record
			mousex,mousey,mousex2,mousey2 : unsigned(9 downto 0);
			btn : std_logic_vector(2 downto 0);
			state_mouse : state_mouse_t;
			minx,miny,maxx,maxy : unsigned(9 downto 0);
			newdata : std_logic;
		end record;
	signal r2reg,r2next : regs2_t;
begin
	floattofixed_unit : entity floattofixed port map(
		a=>fltf_a,operation_nd=>fltf_operation_nd,operation_rfd=>fltf_operation_rfd,
		clk=>clk,sclr=>reset,ce=>'1',result=>fltf_result,rdy=>fltf_rdy
	);
	
	mult_unit : entity mult port map(
		a=>mult_a,b=>mult_b,operation_nd=>mult_operation_nd,
		operation_rfd=>mult_operation_rfd,clk=>clk,sclr=>reset,ce=>'1',
		result=>mult_result,rdy=>mult_rdy
	);
	
	div_unit : entity div port map(
		a=>div_a,b=>div_b,operation_nd=>div_operation_nd,
		operation_rfd=>div_operation_rfd,clk=>clk,sclr=>reset,ce=>'1',
		result=>div_result,rdy=>div_rdy
	);
	
	addsub_unit : entity addsub port map(
		a=>addsub_a,b=>addsub_b,operation=>addsub_operation,operation_nd=>addsub_operation_nd,
		operation_rfd=>addsub_operation_rfd,clk=>clk,sclr=>reset,ce=>'1',
		result=>addsub_result,rdy=>addsub_rdy
	);
	
	fixedtofloat_unit : entity fixedtofloat port map(
		a=>ftf_a,operation_nd=>ftf_operation_nd,operation_rfd=>ftf_operation_rfd,
		clk=>clk,sclr=>reset,ce=>'1',result=>ftf_result,rdy=>ftf_rdy
	);
	
	vga_pixel_unit : entity vga_pixel_gen port map(
			clk=>clk,reset=>reset,pixel_tick=>pixel_tick,rgb=>vpu_rgb,wr_tick=>wr_tick,
			wr_addr=>wr_addr,wr_data=>wr_data,wr_fifo_full=>wr_fifo_full,ad=>ad,
			we_n=>we_n, oe_n=>oe_n,dio_a=>dio_a,ce_a_n=>ce_a_n,ub_a_n=>ub_a_n,lb_a_n=>lb_a_n,
			dio_b=>dio_b,ce_b_n=>ce_b_n,ub_b_n=>ub_b_n,lb_b_n=>lb_b_n,
			video_on=>video_on
		);

	vga_sync_unit : entity vga_sync port map(
			clk=>clk,reset=>reset,hsync=>hsync,vsync=>vsync,video_on=>video_on, p_tick=>pixel_tick,
			pixel_x=>pixel_x,pixel_y=>pixel_y
		);

	mouse_unit : entity mouse port map(
			clk=>clk,reset=>reset,ps2d=>ps2d,ps2c=>ps2c,xm=>xm,ym=>ym,btnm=>btnm,m_done_tick=>m_done_tick
		);

	process(pixel_x,r2reg.mousex,pixel_y,r2reg.mousey,r2reg.mousex2,r2reg.mousey2,vpu_rgb,video_on)
	begin
		rgb <= "000";
		if video_on='1' then
			if (unsigned(pixel_x) = r2reg.mousex) or 
				(unsigned(pixel_y) = r2reg.mousey) or
				(unsigned(pixel_x) = r2reg.mousex2) or 
				(unsigned(pixel_y) = r2reg.mousey2) then
				rgb <= "111";
			else
				rgb <= vpu_rgb;
			end if;
		end if;
	end process;
	
	process(clk,reset)
	begin
		if reset='1' then
			rreg.lineaddress <= (others=>'0');
			rreg.state <= s3;
			rreg.state_reset <= s10;
			rreg.color <= "000";
		elsif rising_edge(clk) then
			rreg <= rnext;
		end if;
	end process;
	
	process(clk,reset)
	begin
		if reset='1' then
			r2reg.mousex <= to_unsigned(200,10);
			r2reg.mousey <= to_unsigned(200,10);
			r2reg.mousex2 <= to_unsigned(200,10);
			r2reg.mousey2 <= to_unsigned(200,10);
			r2reg.minx <= to_unsigned(0,10);
			r2reg.miny <= to_unsigned(0,10);
			r2reg.maxx <= to_unsigned(640,10);
			r2reg.maxy <= to_unsigned(480,10);
			r2reg.state_mouse <= s10;
			r2reg.newdata <= '0';
		elsif rising_edge(clk) then
			r2reg <= r2next;
		end if;
	end process;
	
	sumx <= r2reg.mousex + unsigned(xm(8)&xm);
	sumy <= r2reg.mousey - unsigned(ym(8)&ym);
	
	process(m_done_tick,r2reg,xm,ym,sumx,sumy,btnm,mouse_reset)
	begin
		r2next <= r2reg;
		
		case r2reg.state_mouse is
			when s10 =>
				if m_done_tick = '1' then
					r2next.mousex <= sumx;
					r2next.mousey <= sumy;
					r2next.mousex2 <= sumx;
					r2next.mousey2 <= sumy;
					if btnm(0) = '1' then
						r2next.state_mouse <= s20;
					end if;
				end if;
			when s20 =>
				if m_done_tick = '1' then
					r2next.mousex <= sumx;
					r2next.mousey <= sumy;
					if btnm(0) = '0' then
						r2next.state_mouse <= s30;
					end if;
				end if;
			when s30 =>
				if r2reg.mousex < r2reg.mousex2 then
					r2next.minx <= r2reg.mousex;
					r2next.maxx <= r2reg.mousex2;
				else
					r2next.minx <= r2reg.mousex2;
					r2next.maxx <= r2reg.mousex;
				end if;
				if r2reg.mousey < r2reg.mousey2 then
					r2next.miny <= r2reg.mousey;
					r2next.maxy <= r2reg.mousey2;
				else
					r2next.miny <= r2reg.mousey2;
					r2next.maxy <= r2reg.mousey;
				end if;
				r2next.newdata <= '1';
				r2next.state_mouse <= s10;
		end case;	
		if mouse_reset = '1' then
			r2next.minx <= to_unsigned(0,10);
			r2next.miny <= to_unsigned(0,10);
			r2next.maxx <= to_unsigned(640,10);
			r2next.maxy <= to_unsigned(480,10);
			r2next.newdata <= '0';
		end if;
	end process;
	
	process(rreg,wr_fifo_full,sw,ftf_operation_rfd,ftf_rdy,ftf_result,addsub_operation_rfd,addsub_rdy,
			addsub_result,mult_operation_rfd,mult_rdy,mult_result,div_operation_rfd,div_rdy,div_result,
			fltf_operation_rfd,fltf_rdy,fltf_result,r2reg.minx,r2reg.maxx,r2reg.miny,r2reg.maxy,r2reg.newdata)
	begin
		mouse_reset <= '0';
		led <= "11100000";
		rnext <= rreg;
		wr_addr <= (others=>'0');
		wr_data <= (others=>'0');
		wr_tick <= '0';
		addsub_a <= (others=>'0');
		addsub_b <= (others=>'0');
		addsub_operation <= "000000";
		addsub_operation_nd <= '0';
		ftf_a <= (others=>'0');
		ftf_operation_nd <= '0';
		mult_a <= (others=>'0');
		mult_b <= (others=>'0');
		mult_operation_nd <= '0';
		div_a <= (others=>'0');
		div_b <= (others=>'0');
		div_operation_nd <= '0';
		fltf_a <= (others=>'0');
		fltf_operation_nd <= '0';

		case rreg.state is
			when s3 =>
				case rreg.state_reset is
					-- calculate xmin and ymin
					when s10 =>
						if ftf_operation_rfd = '1' then
							ftf_a <= std_logic_vector(to_signed(-2,12));
							ftf_operation_nd <= '1';
							rnext.state_reset <= s20;
						end if;
					when s20 =>
						if ftf_rdy = '1' then
							rnext.xmin <= ftf_result;
							rnext.ymin <= ftf_result;
							rnext.state_reset <= s30;
						end if;
					-- calculate xmax and ymax
					when s30 =>
						if ftf_operation_rfd = '1' then
							ftf_a <= std_logic_vector(to_signed(2,12));
							ftf_operation_nd <= '1';
							rnext.state_reset <= s40;
						end if;
					when s40 =>
						if ftf_rdy = '1' then
							rnext.xmax <= ftf_result;
							rnext.ymax <= ftf_result;
							rnext.state_reset <= s50;
						end if;
					-- calculate width
					when s50 =>
						if ftf_operation_rfd = '1' then
							ftf_a <= std_logic_vector(to_signed(640,12));
							ftf_operation_nd <= '1';
							rnext.state_reset <= s60;
						end if;
					when s60 =>
						if ftf_rdy = '1' then
							rnext.width <= ftf_result;
							rnext.state_reset <= s70;
						end if;
					-- calculate height
					when s70 =>
						if ftf_operation_rfd = '1' then
							ftf_a <= std_logic_vector(to_signed(480,12));
							ftf_operation_nd <= '1';
							rnext.state_reset <= s80;
						end if;
					when s80 =>
						if ftf_rdy = '1' then
							rnext.height <= ftf_result;
							rnext.state_reset <= s90;
						end if;
					when s90 =>
						rnext.state <= s4;
				end case;
			when s4 =>
				rnext.address <= rreg.lineaddress;
				rnext.state <= s5;
				rnext.x <= (others=>'0');
			when s5 =>
				if wr_fifo_full='0' then
					wr_addr <= std_logic_vector(rreg.address);
					wr_data <= (others=>'1');
					wr_tick <= '1';
					if rreg.address = (rreg.lineaddress + 79) then
						rnext.address <= rreg.lineaddress;
						rnext.state <= calcxm;
						rnext.i <= to_unsigned(8,4);
					else
						rnext.address <= rreg.address + 1;
					end if;
				end if;
			when calcxm =>
				rnext.r1(11 downto 0) <= std_logic_vector(rreg.x);
				rnext.r2 <= rreg.xmax;
				rnext.r3 <= rreg.xmin;
				rnext.r4 <= rreg.width;
				rnext.return_state <= calcym;
				rnext.state <= sub_linear;
				rnext.state_sub_linear <= s5;
			when calcym =>
				rnext.xmath <= rreg.r1;
				
				rnext.r1(11 downto 0) <= std_logic_vector(rreg.y);
				rnext.r2 <= rreg.ymax;
				rnext.r3 <= rreg.ymin;
				rnext.r4 <= rreg.height;
				rnext.return_state <= initcolcal;
				rnext.state <= sub_linear;
				rnext.state_sub_linear <= s5;
			when initcolcal =>
				rnext.ymath <= rreg.r1;
				
				rnext.state <= colorcalculation;
				rnext.state_colcal <= s10;
			when colorcalculation =>
				case rreg.state_colcal is
					when s10 =>
						rnext.r1 <= std_logic_vector(to_unsigned(1,18));
						rnext.color <= "001";
						rnext.r2 <= (others => '0');
						rnext.r3 <= (others => '0');
						rnext.state_colcal <= s20;
					when s20 =>
						if mult_operation_rfd = '1' then
							mult_a <= rreg.r2;
							mult_b <= rreg.r2;
							mult_operation_nd <= '1';
							rnext.state_colcal <= s30;
						end if;
					when s30 =>
						if mult_rdy = '1' then
							rnext.r4 <= mult_result;
							rnext.state_colcal <= s40;							
						end if;
					when s40 =>
						if mult_operation_rfd = '1' then
							mult_a <= rreg.r3;
							mult_b <= rreg.r3;
							mult_operation_nd <= '1';
							rnext.state_colcal <= s50;
						end if;
					when s50 =>
						if mult_rdy = '1' then
							rnext.r5 <= mult_result;
							rnext.state_colcal <= s60;							
						end if;
					when s60 =>
						if addsub_operation_rfd = '1' then
							addsub_a <= rreg.r4;
							addsub_b <= rreg.r5;
							addsub_operation <= "000001";
							addsub_operation_nd <= '1';
							rnext.state_colcal <= s70;
						end if;
					when s70 =>
						if addsub_rdy = '1' then
							rnext.r4 <= addsub_result;
							rnext.state_colcal <= s80;
						end if;
					when s80 =>
						if addsub_operation_rfd = '1' then
							addsub_a <= rreg.r4;
							addsub_b <= rreg.xmath;
							addsub_operation <= "000000";
							addsub_operation_nd <= '1';
							rnext.state_colcal <= s90;							
						end if;
					when s90 =>
						if addsub_rdy = '1' then
							rnext.r4 <= addsub_result;
							rnext.state_colcal <= s100;
						end if;
					when s100 =>
						if ftf_operation_rfd = '1' then
							ftf_a <= "000000000010";
							ftf_operation_nd <= '1';
							rnext.state_colcal <= s110;
						end if;
					when s110 =>
						if ftf_rdy = '1' then
							rnext.r5 <= ftf_result;
							rnext.state_colcal <= s120;
						end if;
					when s120 =>
						if mult_operation_rfd = '1' then
							mult_a <= rreg.r5;
							mult_b <= rreg.r2;
							mult_operation_nd <= '1';
							rnext.state_colcal <= s130;
						end if;
					when s130 =>
						if mult_rdy = '1' then
							rnext.r5 <= mult_result;
							rnext.state_colcal <= s140;
						end if;
					when s140 =>
						if mult_operation_rfd = '1' then
							mult_a <= rreg.r5;
							mult_b <= rreg.r3;
							mult_operation_nd <= '1';
							rnext.state_colcal <= s150;
						end if;
					when s150 =>
						if mult_rdy = '1' then
							rnext.r5 <= mult_result;
							rnext.state_colcal <= s160;
						end if;
					when s160 =>
						if addsub_operation_rfd = '1' then
							addsub_a <= rreg.r5;
							addsub_b <= rreg.ymath;
							addsub_operation <= "000000";
							addsub_operation_nd <= '1';
							rnext.state_colcal <= s170;
						end if;
					when s170 =>
						if addsub_rdy = '1' then
							rnext.r5 <= addsub_result;
							rnext.state_colcal <= s180;
						end if;
					when s180 =>
						rnext.r2 <= rreg.r4;
						rnext.r3 <= rreg.r5;
						rnext.state_colcal <= s190;
					when s190 =>
						if mult_operation_rfd = '1' then
							mult_a <= rreg.r2;
							mult_b <= rreg.r2;
							mult_operation_nd <= '1';
							rnext.state_colcal <= s200;
						end if;
					when s200 =>
						if mult_rdy = '1' then
							rnext.r4 <= mult_result;
							rnext.state_colcal <= s210;
						end if;
					when s210 =>
						if mult_operation_rfd = '1' then
							mult_a <= rreg.r3;
							mult_b <= rreg.r3;
							mult_operation_nd <= '1';
							rnext.state_colcal <= s220;
						end if;
					when s220 =>
						if mult_rdy = '1' then
							rnext.r5 <= mult_result;
							rnext.state_colcal <= s230;
						end if;
					when s230 =>
						if addsub_operation_rfd = '1' then
							addsub_a <= rreg.r4;
							addsub_b <= rreg.r5;
							addsub_operation <= "000000";
							addsub_operation_nd <= '1';
							rnext.state_colcal <= s240;
						end if;
					when s240 =>
						if addsub_rdy = '1' then
							rnext.r6 <= addsub_result;
							rnext.state_colcal <= s250;
						end if;
					when s250 =>
						if fltf_operation_rfd = '1' then
							fltf_a <= rreg.r6;
							fltf_operation_nd <= '1';
							rnext.state_colcal <= s260;
						end if;
					when s260 =>
						if fltf_rdy = '1' then
							rnext.r6 <= "000000" & fltf_result;
							rnext.state_colcal <= s270;
						end if;
					when s270 =>
						if unsigned(rreg.r6) >= 4 then
							rnext.state <= s10;
						elsif unsigned(rreg.r1) = 50 then
							rnext.color <= "000";
							rnext.state <= s10;
						else
							rnext.r1 <= std_logic_vector(unsigned(rreg.r1) + to_unsigned(1,16));
							if rreg.color = "111" then
								rnext.color <= "001";
							else
								rnext.color <= std_logic_vector(unsigned(rreg.color) + to_unsigned(1,3));
							end if;
							rnext.state_colcal <= s60;
						end if;
				end case;
			when s10 =>
				if rreg.i > 0 then
					rnext.dataword <= rreg.dataword(27 downto 0) & rreg.color & '0';
					--rnext.dataword <= "1000" & sw & '0'& sw & '0'& sw & '0'& sw & '0'& sw & '0'& sw & '0'& sw & '0';
					rnext.state <= calcxm;
					rnext.i <= rreg.i - 1;
					rnext.x <= rreg.x + 1;
				elsif wr_fifo_full='0' then
					rnext.i <= to_unsigned(8,4);
					wr_addr <= std_logic_vector(rreg.address);
					wr_data <= rreg.dataword;
					wr_tick <= '1';
					rnext.n <= (others=>'1');
					rnext.state <= s20;
				end if;
			when s20 =>
				if rreg.n = 0 then
					rnext.state <= s22;
				else
					rnext.n <= rreg.n - 1;
				end if;
				rnext.state <= s22;
			when s22 =>
				if rreg.address = (rreg.lineaddress + 82) and r2reg.newdata = '1' then
					rnext.r1(11 downto 0) <= "00" & std_logic_vector(r2reg.minx);
					rnext.r2 <= rreg.xmax;
					rnext.r3 <= rreg.xmin;
					rnext.r4 <= rreg.width;
					rnext.return_state <= s24;
					rnext.state <= sub_linear;
					rnext.state_sub_linear <= s5;					
				else
					rnext.state <= s30;
				end if;
			when s24 =>
				rnext.r5 <= rreg.r1;
				
				rnext.r1(11 downto 0) <= "00" & std_logic_vector(r2reg.maxx);
				rnext.r2 <= rreg.xmax;
				rnext.r3 <= rreg.xmin;
				rnext.r4 <= rreg.width;
				rnext.return_state <= s26;
				rnext.state <= sub_linear;
				rnext.state_sub_linear <= s5;					
			when s26 =>
				rnext.xmax <= rreg.r1;
				rnext.xmin <= rreg.r5;
				
				rnext.r1(11 downto 0) <= "00" & std_logic_vector(r2reg.miny);
				rnext.r2 <= rreg.ymax;
				rnext.r3 <= rreg.ymin;
				rnext.r4 <= rreg.height;
				rnext.return_state <= s28;
				rnext.state <= sub_linear;
				rnext.state_sub_linear <= s5;					
			when s28 =>
				rnext.r5 <= rreg.r1;
				
				rnext.r1(11 downto 0) <= "00" & std_logic_vector(r2reg.maxy);
				rnext.r2 <= rreg.ymax;
				rnext.r3 <= rreg.ymin;
				rnext.r4 <= rreg.height;
				rnext.return_state <= s29;
				rnext.state <= sub_linear;
				rnext.state_sub_linear <= s5;					
			when s29 =>
				rnext.ymax <= rreg.r1;
				rnext.ymin <= rreg.r5;
				mouse_reset <= '1';
				rnext.state <= s30;
			when s30 =>
				if rreg.address = (rreg.lineaddress + 82) then
					if rreg.lineaddress = 38320 then
						rnext.lineaddress <= (others=>'0');
						rnext.y <= (others=>'0');
					else
						rnext.lineaddress <= rreg.lineaddress + 80;
						rnext.y <= rreg.y + 1;
					end if;
					rnext.state <= s4;
				else
					rnext.address <= rreg.address + 1;
					rnext.state <= calcxm;
				end if;
			when sub_linear =>
				case rreg.state_sub_linear is
					when s5 =>
						if ftf_operation_rfd = '1' then
							ftf_a <= rreg.r1(11 downto 0);
							ftf_operation_nd <= '1';
							rnext.state_sub_linear <= s8;
						end if;
					when s8 =>
						if ftf_rdy = '1' then
							rnext.r1 <= ftf_result;
							rnext.state_sub_linear <= s10;
						end if;
					when s10 =>
						if addsub_operation_rfd = '1' then
							addsub_a <= rreg.r2;
							addsub_b <= rreg.r3;
							addsub_operation <= "000001";
							addsub_operation_nd <= '1';
							rnext.state_sub_linear <= s20;
						end if;
					when s20 =>
						if addsub_rdy = '1' then
							rnext.r2 <= addsub_result;
							rnext.state_sub_linear <= s30;							
						end if;
					when s30 =>
						if mult_operation_rfd = '1' then
							mult_a <= rreg.r1;
							mult_b <= rreg.r2;
							mult_operation_nd <= '1';
							rnext.state_sub_linear <= s40;
						end if;
					when s40 =>
						if mult_rdy = '1' then
							rnext.r1 <= mult_result;
							rnext.state_sub_linear <= s50;
						end if;
					when s50 =>
						if div_operation_rfd = '1' then
							div_a <= rreg.r1;
							div_b <= rreg.r4;
							div_operation_nd <= '1';
							rnext.state_sub_linear <= s60;
						end if;
					when s60 =>
						if div_rdy = '1' then
							rnext.r1 <= div_result;
							rnext.state_sub_linear <= s70;
						end if;
					when s70 =>
						if addsub_operation_rfd = '1' then
							addsub_a <= rreg.r1;
							addsub_b <= rreg.r3;
							addsub_operation_nd <= '1';
							addsub_operation <= "000000";
							rnext.state_sub_linear <= s80;
						end if;
					when s80 =>
						if addsub_rdy = '1' then
							rnext.r1 <= addsub_result;
							rnext.state <= rreg.return_state;
						end if;
				end case;
		end case;
	end process;
end Behavioral;
