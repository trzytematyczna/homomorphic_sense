import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.Properties;

import sun.security.krb5.internal.PAData;
import edu.emory.mathcs.jtransforms.dct.FloatDCT_1D;
import edu.emory.mathcs.jtransforms.dct.FloatDCT_2D;
import ij.ImagePlus;
import ij.gui.GenericDialog;
import ij.io.FileSaver;
import ij.io.OpenDialog;
import ij.measure.ResultsTable;
import ij.plugin.GaussianBlur3D;
import ij.plugin.TextReader;
import ij.plugin.filter.GaussianBlur;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;


public class Main {
	
	public static void main(String[] args) {
		
		float[] a = {1, 2, 3, 4, 5, 6, 7, 8, 9};
//		float[] b = {1, 1, 2, 3, 4, 5, 6, 7, 8};
		float[] b = {5, 1, 2, 6, 4, 5, 6, 3, 8};
		float[] c = {1, 1, 2, 3, 4, 5, 6, 7, 8};
		ImageProcessor A = new FloatProcessor(3, 3, a);
		ImageProcessor B = new FloatProcessor(3, 3, b);
//		ImageProcessor O = RiceHommomorfEst_.filter2b(A, B);
		ImageProcessor K = RiceHommomorfEst_.find(B, 4, RiceHommomorfEst_.GOE);
		ImageProcessor O = RiceHommomorfEst_.applyfromfind(B, K, 0);
		for (int i = 0; i < 6; i++) {
			for (int j = 0; j < 1; j++) {
				System.out.print("[" + K.getPixelValue(j, i) + "]");
			}
			System.out.println();
		}
		System.out.println(O.getHeight());
		System.out.println(O.getWidth());
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				System.out.print("[" + O.getPixelValue(j, i) + "]");
			}
			System.out.println();
		}
		System.out.println();
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				System.out.print("[" + B.getPixelValue(j, i) + "]");
			}
			System.out.println();
		}
		
		System.out.println("START");
		
		RiceHommomorfEst_ rhe = new RiceHommomorfEst_();
		rhe.setup(null, null);
		rhe.run(null);
		
		System.out.println("END");
		
		return;
		
	}
	
}
