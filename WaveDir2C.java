import java.util.*;
import java.io.*;
import javax.sound.sampled.*;



public class WaveDir2C {
	static ArrayList<File> files = new ArrayList<>();

	public static void main(String[] args) {

                         for(File f : new File(args[1]).listFiles()) {
                                files.add(f);
                        }
 
		String ws = args[0];
		
		System.out.println("signed short int _" + ws + "s[][200000] = {");
		//System.out.println("let instruments = [");
		
		for(File f : files) {
                System.out.println("{");
		int j = 0;
		try {
                                AudioInputStream audioIn = AudioSystem.getAudioInputStream(f);
				byte[] bs = new byte[1024];
					boolean ok = false;
                                while(audioIn.available() > 0) l : {
                                        int n = audioIn.read(bs);

					for(int i = 0;i < n / 2;++i) {
					j += 1;
					if(j >= 200000) break l;
						int a = (bs[2 * i + 0] & 0xff) | (bs[i * 2 + 1] << 8);
						if(!ok) ok = Math.abs(a) > 100;
						if(ok)
 						System.out.print((a) + ",");
					}

                                }
				System.out.println();


		} catch(Exception e) {
			e.printStackTrace();
		}

                System.out.println("},");

		}

                System.out.println("};");
                
		System.out.println("signed short int *" + ws + "s = &_" + ws + "s[0][0];");
                System.out.println("long num" + ws.substring(0, 1).toUpperCase() + ws.substring(1) + "s = sizeof(_" + ws + "s) / sizeof(_" + ws + "s[0]);");
                System.out.println("long " + ws + "Length = sizeof(_" + ws + "s[0]) / sizeof(_" + ws + "s[0][0]);");


	}
}



