package homomorphic_sense;

import ij.ImagePlus;
import ij.plugin.TextReader;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;


public class Main {
	
	public static void main(String[] args) {
		
//		System.loadLibrary("opencv_java300");
		
		TextReader textReader = new TextReader();
		ImageProcessor mriIp = textReader.open("res/MR_noisy.csv");
		ImageProcessor snrIp = textReader.open("res/MR_SNR.csv");
		System.out.println(mriIp.getPixelValue(0, 0));
		System.out.println(mriIp.getPixelValue(1, 1));
		mriIp.add(4);
		System.out.println(mriIp.getPixelValue(0, 0));
		System.out.println(mriIp.getPixelValue(1, 1));
		ImagePlus asd = new ImagePlus("ASD", mriIp);
	    float[] jedne = {2,1,1,1,1,1};
	    float[] dwa = {1,1,1,1,1,1};
	    float[][] www = {jedne,dwa};
	    float[][] ppp = new float[7][8];
//	    System.out.println(www[3][0]);
	   ImageProcessor asdasd = new FloatProcessor(ppp);
	   System.out.println(asdasd.getPixelValue(0, 0));
	   System.out.println(asdasd.getHeight());
	   System.out.println(asdasd.getWidth());
//		Mat mri = new Mat(mriIp.getHeight(), mriIp.getWidth(), CvType.CV_64FC1);
//		Mat snr = new Mat(snrIp.getHeight(), snrIp.getWidth(), CvType.CV_64FC1);
//		for (int i = 0; i < mriIp.getHeight(); i++){
//			for (int j = 0; j < mriIp.getWidth(); j++) {
//				mri.put(i, j, mriIp.getPixelValue(i, j));
//				snr.put(i, j, snrIp.getPixelValue(i, j));
//			}
//		}
		
//		Mat ricianMap =  RiceHommomorfEst.rice_hommomorf_est(mri, snr, 4.8, 2, RiceHommomorfEst.RICIAN, 3);
//		
//		Imgcodecs.imwrite("MR_Rician_map.csv", ricianMap);
		
	}
}
