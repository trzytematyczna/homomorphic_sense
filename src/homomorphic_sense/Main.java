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
		FFT ff = new FFT();
		ff.fft(mri, false);
		
//		SimpleMatrix gaussian = RiceHomomorfEst.rice_hommomorf_est(mri, snr, 4.8, 2, RiceHomomorfEst.GAUSSIAN);
		SimpleMatrix rician = RiceHomomorfEst.rice_hommomorf_est(mri, snr, 4.8, 2, RiceHomomorfEst.RICIAN, 3.0);
		

	}
	
}
