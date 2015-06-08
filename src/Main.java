import ij.plugin.TextReader;
import ij.process.ImageProcessor;

public class Main {
	
	public static void main(String[] args) {
				
		TextReader textReader = new TextReader();
		ImageProcessor mriIp = textReader.open("res/MR_noisy.csv");
		ImageProcessor snrIp = textReader.open("res/MR_SNR.csv");
		
		FFT_ fft = new FFT_();
		fft.fft(mriIp, false);
		
		Complex c = null;
		
	}
}
