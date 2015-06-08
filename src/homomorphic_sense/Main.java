package homomorphic_sense_opencv;

import ij.plugin.TextReader;
import ij.process.ImageProcessor;

import org.opencv.core.CvType;
import org.opencv.core.Mat;
import org.opencv.imgcodecs.*;
import org.opencv.videoio.*;

public class Main {
	
	public static void main(String[] args) {
		
		System.loadLibrary("opencv_java300");
		
		TextReader textReader = new TextReader();
		ImageProcessor mriIp = textReader.open("res/MR_noisy.csv");
		ImageProcessor snrIp = textReader.open("res/MR_SNR.csv");
		
		Mat mri = new Mat(mriIp.getHeight(), mriIp.getWidth(), CvType.CV_64FC1);
		Mat snr = new Mat(snrIp.getHeight(), snrIp.getWidth(), CvType.CV_64FC1);
		for (int i = 0; i < mriIp.getHeight(); i++){
			for (int j = 0; j < mriIp.getWidth(); j++) {
				mri.put(i, j, mriIp.getPixelValue(i, j));
				snr.put(i, j, snrIp.getPixelValue(i, j));
			}
		}
		
		Mat ricianMap =  RiceHommomorfEst.rice_hommomorf_est(mri, snr, 4.8, 2, RiceHommomorfEst.RICIAN, 3);
		
		Imgcodecs.imwrite("MR_Rician_map.csv", ricianMap);
		
	}
}
