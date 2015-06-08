import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

import com.sun.corba.se.spi.orbutil.fsm.Guard.Result;

import ij.*;
import ij.gui.NewImage;
import ij.measure.ResultsTable;
import ij.plugin.TextReader;
import ij.plugin.filter.PlugInFilter;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import ij.process.FloatProcessor;
import edu.emory.mathcs.jtransforms.fft.*;

public class RiceHommomorfEst_ implements PlugInFilter {

	static int ex_filter_type = 1;     //# 1 - local mean, 
//            # 2 - expectation-maximization (EM).
	static int ex_window_size =5; //# window size for E{X},
	static int ex_iterations= 10;  //   # number of iterations of the EM algorithm (used only by EM),
	static double lpf_f = 3.4;    //        # sigma for LPF filter,
	static double lpf_f_SNR = 1.2;  //      # sigma for LPF filter; used to smooth sigma(x) in SNR,
	static double lpf_f_Rice= 5.4;       //# sigma for LPF filter; used to smooth Rician corrected noise map,
	static String input_filename = "MR_noisy.csv";                    //# Noisy MR image,
	static String input_filenameSNR = "MR_SNR.csv";  //                  # Noisy MR image,
	static String output_filename_Gaussian = "MR_Gaussian_Map.csv";   //# estimated noise map for Gaussian case,
	static String output_filename_Rician = "MR_Rician_Map.csv"  ;//    # estimated noise map for Rician case,
	
	
	static int LT =1;
	static int	GT= 2;
	static int	EQ =3;
	static int	LOE = 4;
	static int	GOE =5;
	static int	NEQ =6;
	
	@Override
	public void run(ImageProcessor arg0) {
		TextReader textReader = new TextReader();
		ImageProcessor mriIp = textReader.open("res/MR_noisy.csv");
		ImageProcessor snrIp = textReader.open("res/MR_SNR.csv");
		ImageProcessor[] res = rice_hommomorf_est(mriIp, snrIp, RiceHommomorfEst_.lpf_f, ex_filter_type,??, ex_window_size);
		ResultsTable rician =  ResultsTable.createTableFromImage(res[0]);
		ResultsTable gauss =  ResultsTable.createTableFromImage(res[1]);
		rician.saveAs(output_filename_Rician);
		gauss.saveAs(output_filename_Gaussian);
	}

	@Override
	public int setup(String arg0, ImagePlus arg1) {
		// TODO Auto-generated method stub
		
		BufferedReader br = null;
		try {
 
			String sCurrentLine;
			br = new BufferedReader(new FileReader("C:\\testing.txt"));
			while ((sCurrentLine = br.readLine()) != null) {
				System.out.println(sCurrentLine);
			}
		} catch (IOException e) {
			e.printStackTrace();
		} finally {
			try {
				if (br != null)br.close();
			} catch (IOException ex) {
				ex.printStackTrace();
			}
		}
		
		return  DOES_ALL+NO_IMAGE_REQUIRED;
	}
	
	public static ImageProcessor[] rice_hommomorf_est(ImageProcessor In, ImageProcessor SNR, double LPF, int Modo, int noiseType, int winsize) {
		
		ImageProcessor[] result = new FloatProcessor[2];
		
		ImageProcessor MapaR = null;
		ImageProcessor MapaG = null;
		
		int[] Ws = {winsize, winsize};
		
		ImageProcessor [] em_ml = em_ml_rice2D(In, 10, Ws);
		ImageProcessor M2 = em_ml[0];
		ImageProcessor Sigma_n = em_ml[1];
		ImageProcessor M1;
		
		
		Sigma_n = lpf(Sigma_n,1.2);
		
		M1 = filter2b(createImage(5, 5, 0.25), In);
		
		if(SNR.getHeight() == 1 && SNR.getHeight() == 0){
			SNR = divide(M2, Sigma_n);
		}
		
		ImageProcessor LocalMean;
		ImageProcessor Rn ;
		ImageProcessor LPF1;
		ImageProcessor LPF2;
		ImageProcessor Fc1;

		double psi = -0.5772;
		double exp_psi_div2 = 1.3346;

		//RICIAN!!
//		if(noiseType == RICIAN){
			if(Modo == 1){
				LocalMean = M1;
			}
			else if (Modo == 2){
				LocalMean = M2;
			}
			else{
				LocalMean = createImage(In.getWidth(), In.getHeight(), 0);
			}
			
			Rn = absdiff(In, LocalMean);
					
			ImageProcessor lRn = add(multiply(Rn, compare(Rn, 0, NEQ)), multiply(compare(Rn, 0, EQ),0.001));
			lRn.log();
			LPF2 = lpf(lRn,LPF);
			Fc1 = correct_rice_gauss(SNR);
			LPF1 = substract(LPF2,Fc1);
			LPF1 = lpf(LPF1, LPF+2,2);
			LPF1.exp();
			ImageProcessor Mapa1 = LPF1;

			MapaR = multiply(Mapa1, 2);
			MapaR.sqrt();
			MapaR = multiply(MapaR, exp_psi_div2);
//		}
		
			//GAUSSIAN!!
//		else if (noiseType == GAUSSIAN){
			Rn = absdiff(In, M1);
			ImageProcessor lRn = add(multiply(Rn, compare(Rn, 0, NEQ)), multiply(compare(Rn, 0, EQ),0.001));
			lRn.log();
			LPF2 = lpf(lRn,LPF);
			LPF2.exp();
			ImageProcessor Mapa2 = LPF2;
			MapaG = multiply(Mapa2, 2);
			MapaG.sqrt();
			MapaG = multiply(MapaG, exp_psi_div2);
			
//		}
		result[0] = MapaR;
		result[1] = MapaG;
		
		return result;
	}

	private static ImageProcessor correct_rice_gauss(ImageProcessor SNR) {
		//Fc=Coefs(1)+Coefs(2).*a1+Coefs(3).*a1.^2+Coefs(4).*a1.^3+Coefs(5).*a1.^4+Coefs(6).*a1.^5+Coefs(7).*a1.^6+Coefs(8).*a1.^7+Coefs(9).*a1.^8;
		//Fc=Fc.*(a1<=7);       

		double[] Coefs = {-0.289549906258443,	-0.0388922575606330,	0.409867108141953,	-0.355237628488567,
				0.149328280945610,	-0.0357861117942093,	0.00497952893859122,	-0.000374756374477592,	1.18020229140092e-05};

		ImageProcessor Fc = add(multiply(pow(SNR, 8), Coefs[8]), multiply(pow(SNR, 7), Coefs[7]));
		Fc = add(Fc, add(multiply(pow(SNR, 6), Coefs[6]), multiply(pow(SNR, 5), Coefs[5])));
		Fc = add(Fc, add(multiply(pow(SNR, 4), Coefs[4]), multiply(pow(SNR, 3), Coefs[3])));
		Fc = add(Fc, add(multiply(pow(SNR, 2), Coefs[2]), multiply(pow(SNR, 1), Coefs[1])));
		Fc = add(Fc, Coefs[0]);
		
		Fc = multiply(Fc, compare(SNR, 7, LOE));

		return Fc;
	}

	private static ImageProcessor filter2b(ImageProcessor h, ImageProcessor I) {
		int Mx = h.getHeight();
		int My = h.getWidth();
//		System.out.println(Mx);
//		System.out.println(My);
		
		if(Mx % 2 == 0 || My % 2 == 0){
			System.err.println("filter2b size of h must be odd");
			return null;
		}
		
		int Nx = (Mx - 1) / 2;
		int Ny = (My - 1) / 2;
		ImageProcessor It = padarray(I, Nx, Ny);
		ImageProcessor I2 = filter2(h, It);
//		I_out=I2((Nx+1):end-Nx,(Ny+1):end-Ny);
		ImageProcessor I_out = submatrix(I2, Ny + 1, I2.getWidth() - Ny,  Nx + 1, I2.getHeight() - Nx);
		
		return I_out;
		
	}

	private static ImageProcessor submatrix(ImageProcessor mat, int beginW, int endW,int beginH, int endH) {
		float[][] outf = new float[endW-beginW + 1][endH-beginH + 1];
		for (int i = 0; i < endW-beginW + 1; i++) {
			for (int j = 0; j < endH-beginH + 1; j++) {
				outf[i][j]= mat.getPixelValue(i + beginW,j + beginH);
			}
		}
		ImageProcessor out = new FloatProcessor(outf);
		return out;
	}

	private static ImageProcessor padarray(ImageProcessor in, int nx, int ny) { //chyba done
		int rows = in.getHeight() + nx * in.getHeight();
		int cols = in.getWidth() + ny * in.getWidth();
		float[][] outf = new float[rows][cols];
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < cols; j++) {
				float value = in.get(i / (1 + nx), j / (1 + ny));
				outf[i][j]= value;
			}
		}
		ImageProcessor output = new FloatProcessor(outf);
		return output;
	}

	private static ImageProcessor lpf(ImageProcessor sigma_n, double d) {
		// TODO Auto-generated method stub
		return null;
	}

	private static ImageProcessor lpf(ImageProcessor I, double sigma, int MODO) {
		//TODO

		if (MODO == 1) {
			int Mx = I.getHeight();
			int My = I.getWidth();
			ImageProcessor h = Imgproc.getGaussianKernel(I.getHeight() , sigma, CvType.CV_64F);
			h = divide(h, max(h));
			if (Mx == 1 ||My == 1) {
				ImageProcessor lRnF = new ImageProcessor(I.size(), I.type());
				Core.dft(I, lRnF);
				fftshift(lRnF);
				ImageProcessor lRnF2 = multiply(lRnF, h);
				fftshift(lRnF2);
				ImageProcessor If = new ImageProcessor(I.size(), I.type());
				Core.idft(lRnF2, If);
				return If;
			}
			else {
				ImageProcessor lRnF = new ImageProcessor(I.size(), I.type());
				Core.dft(I, lRnF);
				fftshift(lRnF);
				ImageProcessor lRnF2 = multiply(lRnF, h);
				fftshift(lRnF2);
				ImageProcessor If = new ImageProcessor(I.size(), I.type());
				Core.idft(lRnF2, If);
				return If;
			}
		} else if (MODO == 2) {
			int Mx = I.getHeight();
			int My = I.getWidth();
			ImageProcessor h = Imgproc.getGaussianKernel(I.getHeight() * 2 , sigma * 2, CvType.CV_64F);
			h = divide(h, max(h));
			h = h.submat(new Range(Mx + 1, h.getHeight()), new Range(My + 1, h.getWidth()));
			if (Mx == 1 ||My == 1) {
				ImageProcessor lRnF = new ImageProcessor(I.size(), I.type());
				Core.dct(I, lRnF);
				ImageProcessor lRnF2 = multiply(lRnF, h);
				ImageProcessor If = new ImageProcessor(I.size(), I.type());
				Core.idct(lRnF2, If);
				return If;
			}
			else {
				ImageProcessor lRnF = new ImageProcessor(I.size(), I.type());
				Core.dct(I, lRnF);
				ImageProcessor lRnF2 = multiply(lRnF, h);
				ImageProcessor If = new ImageProcessor(I.size(), I.type());
				Core.idct(lRnF2, If);
				return If;
			}
		}
		return null;
	}
	
	private static ImageProcessor approxI1_I0(ImageProcessor z) {
		//TODO: finish

//		cont=(z<1.5);
		ImageProcessor cont = compare(z, 1.5, RiceHommomorfEst_.LOE);
//		z8=8.*z;
		ImageProcessor z8 = multiply(z, 8);

//		Mn=1-3./z8-15./2./(z8).^2-(3*5*21)./6./(z8).^3;
		ImageProcessor Mn = substract(substract(substract(1, divide(3, z8)), divide(15/2, pow(z8, 2))), divide((3*5*21)/6, pow(z8, 3)));
		
//		Md=1+1./z8+9./2./(z8).^2+(25*9)./6./(z8).^3;
		ImageProcessor Md = add(add(add(divide(1, z8), 1), divide(9/2, z8)), divide((25*9)/6, pow(z8, 3)));
		
//		M=Mn./Md;
		ImageProcessor M = divide(Mn, Md);
//
//		%K=find(isnan(M));
//		%M(K)=besseli(1,z(K))./besseli(0,z(K));
//
//		if (sum(cont)>1)
//		K=find(z<1.5);
		ImageProcessor K = null;
		if(sum(cont)>1){
			K = find(z, 1.5, LT);			
			M = applyfromfind(M, K, divide(besseli(1, valuesfromfind(z,K)), besseli(0, valuesfromfind(z,K))));
		}
		K = find(z,0, EQ);
		
		M = applyfromfind(M, K,0.0);
		
//		M(K)=besseli(1,z(K))./besseli(0,z(K));
//		end
//		K=find(z==0);
//		M(K)=0;
		return M;
	}

	private static ImageProcessor besseli(int besval, ImageProcessor mat) {
		float[][] outf = new float[mat.getWidth()][mat.getHeight()]; 
		if(besval == 0){
			for(int i=0; i<mat.getWidth();i++){
				for(int j=0; j<mat.getHeight();j++){
					outf[i][j] = (float) Besseli.bessj0(mat.getPixelValue(i, j));
				}
			}
		}
		if(besval == 1){
			for(int i=0; i<mat.getWidth();i++){
				for(int j=0; j<mat.getHeight();j++){
					outf[i][j] =  (float) Besseli.bessj1(mat.getPixelValue(i, j));
				}
			}
		}
		return new FloatProcessor(outf);
	}

	private static ImageProcessor valuesfromfind(ImageProcessor z, ImageProcessor k) {
		float[][] outf = new float[1][k.getHeight()];
		int linear_index = 1;
		int iterator = 0;
		
		for(int i=0; i<z.getWidth();i++){
			for(int j=0; j<z.getHeight();j++){
				if(linear_index == k.getPixelValue(1, iterator)){
					outf[1][iterator] = z.getPixelValue(i, j);					
					iterator++;
				}
				linear_index++;
			}
		}		
		return new FloatProcessor(outf);
	}

	private static ImageProcessor applyfromfind(ImageProcessor mat,ImageProcessor k, ImageProcessor z) {
		int linear_index = 1;
		int iterator = 0;
		
		float[][] outf= new float[mat.getWidth()][mat.getHeight()];
		for(int i=0; i<mat.getWidth();i++){
			for(int j=0; j<mat.getHeight();j++){
				if(linear_index == k.getPixelValue(1, iterator)){
					outf[i][j] = z.getPixelValue(1, iterator);					
					iterator++;
				}
				else{
					outf[i][j] = mat.getPixelValue(i, j);
				}
				linear_index++;
			}
		}
		return new FloatProcessor(outf);
	}
	
	private static ImageProcessor applyfromfind(ImageProcessor mat,ImageProcessor k, double val) {
		int linear_index = 1;
		int iterator = 0;
		
		float[][] outf= new float[mat.getWidth()][mat.getHeight()];
		for(int i=0; i<mat.getWidth();i++){
			for(int j=0; j<mat.getHeight();j++){
				if(linear_index == k.getPixelValue(1, iterator)){
					iterator++;
					outf[i][j] = (float) val;					
				}
				else{
					outf[i][j] = mat.getPixelValue(i, j);
				}
				linear_index++;
			}
		}
		return new FloatProcessor(outf);
	}

	private static ImageProcessor find(ImageProcessor mat, double value, int cmpop) {
		int linear_index = 1;
		int iterator = 0;
		int[] temp = new int[mat.getWidth()*mat.getHeight()];
		for(int i=0; i<mat.getWidth();i++){
			for(int j=0; j<mat.getHeight();j++){
				switch (cmpop) {
		            case 1:  if (mat.getPixelValue(i, j) < value){//		int LT =1;
		            			temp[iterator++]=linear_index;
		            		 }
		            		 break;
		            case 2:  if (mat.getPixelValue(i, j) > value){//		int	GT= 2;
		            			temp[iterator++]=linear_index;
		            		 }
		            		 break;
		            case 3:  if (mat.getPixelValue(i, j) == value){//		int	EQ =3;
		            			temp[iterator++]=linear_index;
		            		 }
		            		 break;
		            case 4:  if (mat.getPixelValue(i, j) <= value){//		int	LOE = 4;
		            			temp[iterator++]=linear_index;
		            		 }
		            		 break;
		            case 5:  if (mat.getPixelValue(i, j) >= value){//		int	GOE =5;
		            			temp[iterator++]=linear_index;
		            		 }
		            		 break;
					case 6:  if (mat.getPixelValue(i, j) != value){//		int	NEQ =6;
								temp[iterator++]=linear_index;
				    		 }
				    		 break;
				}
				
				linear_index++;		 
		   }
		}
		float[][] outf = new float[1][iterator];
		for(int i = 0; i<iterator; i++){
			outf[1][i] = temp[i];
		}
		
		return new FloatProcessor(outf);
	}

	private static double sum(ImageProcessor mat1) {
		double res = 0.0;
		for(int i=0; i<mat1.getWidth();i++){
			for(int j=0; j<mat1.getHeight();j++){
				res += mat1.getPixelValue(i, j);
			}			
		}
		return res;
	}

	public static ImageProcessor[] em_ml_rice2D(ImageProcessor In, int N, int[] Ws) {
		
		int Mx = In.getHeight();
		int My = In.getWidth();
		
		float prod = Ws[0] * Ws[1];
		ImageProcessor Mask = divide(createImage(Ws[1],Ws[0], 1.0), prod);
		
//		a_k=sqrt(sqrt(max(2.*filter2B(Mask,In.^2).^2-filter2B(Mask,In.^4),0)));
		ImageProcessor a_k = sqrt(max(substract(multiply(pow(filter2b(Mask, pow(In, 2)), 2), 2), filter2b(Mask, pow(In, 4))), 0));
		
//		sigma_k2=0.5.*max(filter2B(Mask,In.^2)-a_k.^2,0.01);
		ImageProcessor sigma_k2 = multiply(max(substract(filter2b(Mask, pow(In, 2)), pow(a_k, 2)), 0.01), 0.5);
		
		for(int i=1; i<N; i++){
			
//			a_k=max(filter2B(Mask,approxI1_I0(a_k.*In./sigma_k2).*In),0);
			a_k = max(filter2b(Mask, multiply(approxI1_I0(divide(multiply(a_k, In), sigma_k2)), In)), 0);
			
//			sigma_k2=max(0.5.*filter2B(Mask,abs(In).^2)-a_k.^2./2,0.01);
			sigma_k2 = max(substract(multiply(filter2b(Mask, pow(abs(In), 2)), 0.5), divide(pow(a_k, 2), 2)), 0.01);
			
		}
			
		
		ImageProcessor[] result = new ImageProcessor[2];
		result[0] = a_k;
		result[1] = sqrt(sigma_k2);
		
		return result;
	}
	
	private static ImageProcessor fftshift(ImageProcessor mat) {
						
		float[][] shift = new float[mat.getWidth()][mat.getHeight()];
		float[][] pixelsCopy = mat.getFloatArray();
		
		for (int i = 0; i < mat.getWidth() / 2; i++) {
			for (int j = 0; j < mat.getHeight() / 2; j++) {
				shift[i][j] = pixelsCopy[i + mat.getWidth() / 2][j + mat.getHeight() / 2];
				shift[i + mat.getWidth() / 2][j] = pixelsCopy[i][j + mat.getHeight() / 2];
				shift[i][j + mat.getHeight() / 2] = pixelsCopy[i + mat.getWidth() / 2][j];
				shift[i + mat.getWidth() / 2][j + mat.getHeight() / 2] = pixelsCopy[i][j];
			}
		}
		
		ImageProcessor ipshift = new FloatProcessor(shift);
		return ipshift;
		
	}
	
	private static ImageProcessor fft(ImageProcessor ip) {
		float[] pixelsAndZeros = new float[ip.getHeight() * ip.getWidth() * 2];
		float[] pixels = (float[])ip.getPixels();
		for (int i = 0; i < ip.getHeight() * ip.getWidth(); i++) {
			pixelsAndZeros[i] = pixels[i];
		}
		for (int i = ip.getHeight() * ip.getWidth(); i < ip.getHeight() * ip.getWidth() * 2; i++) {
			pixelsAndZeros[i] = 0;
		}
		FloatFFT_2D fft2d = new FloatFFT_2D(ip.getHeight(), ip.getWidth());
		fft2d.realForwardFull(pixelsAndZeros);
		ImageProcessor fftout = new FloatProcessor(ip.getWidth() * 2, ip.getHeight(), pixelsAndZeros);
		return fftout;
	}
	
	private static ImageProcessor ifft(ImageProcessor ip) {
		float[] pixelsCopy = (float[])ip.getPixelsCopy();
		FloatFFT_2D fft2d = new FloatFFT_2D(ip.getHeight(), ip.getWidth() / 2);
		fft2d.complexInverse(pixelsCopy, true);
		float[] pixels = new float[ip.getHeight() * ip.getWidth() / 2];
		for (int i = 0; i < ip.getHeight() * ip.getWidth() / 2; i++) {
			pixels[i] = pixelsCopy[i * 2];
		}
		ImageProcessor ifftout = new FloatProcessor(ip.getWidth() / 2, ip.getHeight(), pixels);
		return ifftout;
	}
	
	private static ImageProcessor createImage(int width, int height, double fill) {
		float[][] data = new float[width][height];
		for (int i = 0; i < width; i++) {
			for (int j = 0; j < height; j++) {
				data[i][j] = (float) fill;
			}
		}
		ImageProcessor ip = new FloatProcessor(data);
		return ip;
	}
	
	private static ImageProcessor absdiff(ImageProcessor mat1, ImageProcessor mat2) {
		return abs(substract(mat1, mat2));
	}
	

	private static ImageProcessor filter2(ImageProcessor mask, ImageProcessor mat) {
		ImageProcessor kernel= flip(mask);
		ImageProcessor out = new FloatProcessor(Convolution.convolute(mat, kernel));
		return out;
	}
	
	
	private static ImageProcessor flip(ImageProcessor mask) {
		// TODO Auto-generated method stub
		return null;
	}

	private static ImageProcessor add(ImageProcessor mat1, ImageProcessor mat2) { //done
		float[][] outf = new float[mat1.getWidth()][mat1.getHeight()];
		for(int i=0; i<mat1.getWidth();i++){
			for(int j=0; j<mat1.getHeight();j++){
				outf[i][j] = mat1.getPixelValue(i, j)+mat2.getPixelValue(i, j);
			}			
		}
		ImageProcessor out = new FloatProcessor(outf);
		return out;
	}
	
	private static ImageProcessor add(ImageProcessor mat, double value) {//done
		mat.add(value);
		return mat;
	}
	
	private static ImageProcessor abs(ImageProcessor mat) {//done
		mat.abs();
		return mat;
	}
	
	private static ImageProcessor substract(ImageProcessor mat1, ImageProcessor mat2) { //done
		float[][] outf = new float[mat1.getWidth()][mat1.getHeight()];
		for(int i=0; i<mat1.getWidth();i++){
			for(int j=0; j<mat1.getHeight();j++){
				outf[i][j] = mat1.getPixelValue(i, j)-mat2.getPixelValue(i, j);
			}			
		}
		ImageProcessor out = new FloatProcessor(outf);
		return out;
	}
	
	private static ImageProcessor substract(double value, ImageProcessor mat) {//done
		ImageProcessor newip = substract(createImage(mat.getWidth(),mat.getHeight(), value), mat);
		return newip ;
	}
	
	private static ImageProcessor sqrt(ImageProcessor mat) {//done
		mat.sqrt();
		return mat;
	}
	
	private static ImageProcessor pow(ImageProcessor mat, double power) {//done
		for(int i=0; i<power; i++){
			mat = multiply(mat, mat);			
		}
		return mat;
	}
	
	private static double max(ImageProcessor mat) { //done
		float max = Integer.MIN_VALUE;
		for(int i=0; i<mat.getWidth();i++){
			for(int j=0; j<mat.getHeight();j++){
				float val =mat.getPixelValue(i, j);
				if(max < val){
					max = val;
				}
			}			
		}
		return max;
	}
	
	private static ImageProcessor max(ImageProcessor mat, double d) { //chyba done
		float[][] outf = new float[mat.getWidth()][mat.getHeight()];
		for(int i=0; i<mat.getWidth();i++){
			for(int j=0; j<mat.getHeight();j++){
				float matval =mat.getPixelValue(i, j);	
				if(d > matval){ //wstawiam wartosc value
					outf[i][j] = (float) d;
				}
				else{ //wstawiam ten sam element co byl
					outf[i][j] = matval;
				}
			}
		}
		ImageProcessor newip = new FloatProcessor(outf);
		return newip;
	}
	
	private static double min(ImageProcessor mat) { //done
		float min = Integer.MAX_VALUE;
		for(int i=0; i<mat.getWidth();i++){
			for(int j=0; j<mat.getHeight();j++){
				float val =mat.getPixelValue(i, j);
				if(min > val){
					min = val;
				}
			}			
		}
		return min;
	}
	
	private static ImageProcessor divide(ImageProcessor mat, double value) {	//done
		mat.multiply(1/value);
		return mat;
	}
	
	private static ImageProcessor divide(double value, ImageProcessor mat) {	//done
		ImageProcessor newip = divide(createImage(mat.getWidth(),mat.getHeight(), value), mat);
		return newip ;
	}
	
	private static ImageProcessor divide(ImageProcessor mat1, ImageProcessor mat2) {	//done
		float[][] outf = new float[mat1.getWidth()][mat1.getHeight()];
		for(int i=0; i<mat1.getWidth();i++){
			for(int j=0; j<mat1.getHeight();j++){
				outf[i][j] = mat1.getPixelValue(i, j)/mat2.getPixelValue(i, j);
			}			
		}
		ImageProcessor out = new FloatProcessor(outf);
		return out;
	}
	
	private static ImageProcessor compare(ImageProcessor mat, double value, int cmpop) { //done
		float[][] outf = new float[mat.getWidth()][mat.getHeight()];
		for(int i=0; i<mat.getWidth();i++){
			for(int j=0; j<mat.getHeight();j++){
				switch (cmpop) {
		            case 1:  if (mat.getPixelValue(i, j) < value){//		int LT =1;
		            			outf[i][j] = 1;
		            		 }
		            		 else {
		         				outf[i][j] = 0;
		            		 }
		            		 break;
		            case 2:  if (mat.getPixelValue(i, j) > value){//		int	GT= 2;
		            			outf[i][j] = 1;
		            		 }
		            		 else {
		         				outf[i][j] = 0;
		            		 }
		            		 break;
		            case 3:  if (mat.getPixelValue(i, j) == value){//		int	EQ =3;
		            			outf[i][j] = 1;
		            		 }
		            		 else {
		         				outf[i][j] = 0;
		            		 }
		            		 break;
		            case 4:  if (mat.getPixelValue(i, j) <= value){//		int	LOE = 4;
		            			outf[i][j] = 1;
		            		 }
		            		 else {
		         				outf[i][j] = 0;
		            		 }
		            		 break;
		            case 5:  if (mat.getPixelValue(i, j) >= value){//		int	GOE =5;
		            			outf[i][j] = 1;
		            		 }
		            		 else {
		         				outf[i][j] = 0;
		            		 }
		            		 break;
					case 6:  if (mat.getPixelValue(i, j) != value){//		int	NEQ =6;
				    			outf[i][j] = 1;
				    		 }
				    		 else {
				 				outf[i][j] = 0;
				    		 }
				    		 break;
					}
				    		 
		        }
		}			
		ImageProcessor out = new FloatProcessor(outf);
		return out;
	}
	
	private static ImageProcessor multiply(ImageProcessor mat, double value) {	//done
		mat.multiply(value);
		return mat;
	}
	private static ImageProcessor multiply(ImageProcessor mat1, ImageProcessor mat2) {	//done
		float[][] outf = new float[mat1.getWidth()][mat1.getHeight()];
		for(int i=0; i<mat1.getWidth();i++){
			for(int j=0; j<mat1.getHeight();j++){
				outf[i][j] = mat1.getPixelValue(i, j)*mat2.getPixelValue(i, j);
			}			
		}
		ImageProcessor out = new FloatProcessor(outf);
		return out;
	}
	
}
