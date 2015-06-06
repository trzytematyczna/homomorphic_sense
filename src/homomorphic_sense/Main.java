package homomorphic_sense;

import org.ejml.equation.Equation;
import org.ejml.simple.SimpleMatrix;

import ij.*;
import ij.io.Opener;
import ij.plugin.TextReader;
import ij.process.ImageProcessor;

public class Main {

	public static void main(String[] args) {
		
		TextReader textReader = new TextReader();
		
		ImageProcessor mriIp = textReader.open("res/MR_noisy.csv");
		ImageProcessor snrIp = textReader.open("res/MR_SNR.csv");
		
		double[][] mriD = new double[mriIp.getHeight()][mriIp.getWidth()];
		double[][] snrD = new double[snrIp.getHeight()][snrIp.getWidth()];
		for (int i = 0; i < mriIp.getHeight(); i++){
			for (int j = 0; j < mriIp.getWidth(); j++) {
				mriD[i][j] = mriIp.getPixelValue(i, j);
				snrD[i][j] = snrIp.getPixelValue(i, j);
			}
		}
		
		SimpleMatrix mri = new SimpleMatrix(mriD);
		SimpleMatrix snr = new SimpleMatrix(snrD);

//		double[] a1 = {1, 2};
//		double[] a2 = {3, 4};
//		double[][] d1 = {a1,a2};
//		double[][] d2 = {a1,a1};
//		SimpleMatrix x = new SimpleMatrix(d1);
//		SimpleMatrix y = new SimpleMatrix(d2);
//		System.out.println(x.get(1, 0));
//		System.out.println(y.get(1,0));
//		System.out.println(x.elementMult(y));
//		System.out.println(x.mult(y));
		
//		SimpleMatrix gaussian = RiceHomomorfEst.rice_hommomorf_est(mri, snr, 4.8, 2, RiceHomomorfEst.GAUSSIAN);
//		SimpleMatrix rician = RiceHomomorfEst.rice_hommomorf_est(mri, snr, 4.8, 2, RiceHomomorfEst.RICIAN);
		

	}
	
}
