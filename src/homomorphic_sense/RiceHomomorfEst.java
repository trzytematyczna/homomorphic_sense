package homomorphic_sense;

import java.util.LinkedList;

import org.ejml.equation.Equation;
import org.ejml.equation.ManagerTempVariables;
import org.ejml.equation.Operation;
import org.ejml.equation.Operation.Info;
import org.ejml.equation.VariableDouble;
import org.ejml.equation.VariableInteger;
import org.ejml.simple.SimpleMatrix;

import ij.process.ImageProcessor;

public class RiceHomomorfEst {

	static int GAUSSIAN = 1;
	static int RICIAN = 2;
	static double[] coeffs = {-0.289549906258443,	-0.0388922575606330,	0.409867108141953,	
		-0.355237628488567,	0.149328280945610,	-0.0357861117942093,	
		0.00497952893859122,	-0.000374756374477592,	1.18020229140092e-05}; 
	
	public static SimpleMatrix rice_hommomorf_est(SimpleMatrix In, SimpleMatrix SNR, double LPF, int Modo, int noiseType) {
		
		double [] Ws = {3,3};
		SimpleMatrix [] em_ml = em_ml_rice2D(In, 10, Ws);
		SimpleMatrix M2 = em_ml[0];
		SimpleMatrix Sigma_n = em_ml[1];
		SimpleMatrix M1;
		
		Sigma_n = lpf(Sigma_n,1.2);
		
		M1 = filter2b(makeOnes(5,5), In);
		
		if(SNR.numRows() == 1 && SNR.numRows() == 0){
			SNR = M2.elementDiv(Sigma_n);
		}
		
		SimpleMatrix LocalMean;
		SimpleMatrix Rn;
		SimpleMatrix LPF1;
		SimpleMatrix LPF2;
		SimpleMatrix Fc1;

		if(noiseType == RICIAN){
			if(Modo == 1){
				LocalMean = M1;
			}
			else if (Modo == 2){
				LocalMean = M2;
			}
			else{
				double[][] zero = new double[1][1];
				zero[0][0] = 0;
				LocalMean = new SimpleMatrix(zero);
			}
			Rn = In.minus(LocalMean);
			SimpleMatrix lRn = new SimpleMatrix(Rn.numRows(), Rn.numCols());
			
			for (int i=0; i<Rn.numRows(); i++){
				for (int j=0; j<Rn.numCols(); j++){
					Rn.set(i, j, Math.abs(Rn.get(i, j)));
				}
			}
			
			for (int i=0; i<Rn.numRows(); i++){
				for (int j=0; j<Rn.numCols(); j++){
					if(Rn.get(i,j) != 0 ){
						lRn.set(i,j, Rn.get(i,j)); 						
					}
					else{
						lRn.set(i,j,0.001);
					}
				}
			}
			
			lRn.elementLog();
			LPF2 = lpf(lRn,LPF);
			Fc1 = correct_rice_gauss(SNR);
			LPF1 = LPF2.minus(Fc1);
			LPF1 = lpf(LPF1, LPF+2,2);
			SimpleMatrix Mapa1 = LPF1.elementExp();
			
			double psi = -0.5772;
			for (int i=0; i<Mapa1.numRows(); i++){
				for (int j=0; j<Mapa1.numCols(); j++){
					double val = Mapa1.get(i, j)*2/Math.sqrt(2)*Math.exp(-psi/2);
					Mapa1.set(i, j, val);
				}
			}
			
			return Mapa1;
		}
		else if (noiseType == GAUSSIAN){
			
		}
		
		return null;
		
	}
	
	private static SimpleMatrix makeOnes(int H, int W) {
		SimpleMatrix res = new SimpleMatrix(W, H);
		for (int i=0; i<H; i++){
			for (int j=0; j<W; j++){
				res.set(i,j, 1);
			}
		}
		return res;
	}

	private static SimpleMatrix lpf(SimpleMatrix lPF1, double d, int i) {
		// TODO Auto-generated method stub
		return null;
	}

	private static SimpleMatrix correct_rice_gauss(SimpleMatrix SNR) {

		//Fc=Coefs(1)+Coefs(2).*a1+Coefs(3).*a1.^2+Coefs(4).*a1.^3+Coefs(5).*a1.^4+Coefs(6).*a1.^5+Coefs(7).*a1.^6+Coefs(8).*a1.^7+Coefs(9).*a1.^8;
		//
		//Fc=Fc.*(a1<=7);       

		SimpleMatrix Fc = new SimpleMatrix(SNR.numRows(), SNR.numCols());
		Fc = multipleM(SNR, coeffs[1]);
		Fc.plus(coeffs[0]).plus(multipleM(SNR.elementPower(2), coeffs[2])).plus(multipleM(SNR.elementPower(3), coeffs[3])).plus(multipleM(SNR.elementPower(4), coeffs[4])).plus(multipleM(SNR.elementPower(5), coeffs[5])).plus(multipleM(SNR.elementPower(6), coeffs[6])).plus(multipleM(SNR.elementPower(7), coeffs[7])).plus(multipleM(SNR.elementPower(8), coeffs[8]));
		
		for (int i=0; i<Fc.numRows(); i++){
			for (int j=0; j<Fc.numCols(); j++){
				double snr_val = SNR.get(i, j);
				if(snr_val <=7){
					Fc.set(i, j, Fc.get(i, j)*snr_val);
				}
			}
		}
		
		return Fc;
	}

	private static SimpleMatrix filter2b(SimpleMatrix ones, SimpleMatrix in) {
//		[Mx, My]=size(h);
//		if (rem(Mx,2)==0)||(rem(My,2)==0)
//		        error('h size must be odd');
//		end
//
//		Nx=(Mx-1)/2;
//		Ny=(My-1)/2;
//		It=padarray(I, [Nx,Ny], 'replicate');
//		%It=im_expand(I,Nx,Ny);
//
//		I2=filter2(h,It);
//		I_out=I2((Nx+1):end-Nx,(Ny+1):end-Ny);		
		if(ones.numRows() % 2 == 0 || ones.numCols() % 2 == 0){
			System.err.println("filter2b size of h must be odd");
			return null;
		}
		
		int Nx = (ones.numRows()-1)/2;
		int Ny = (ones.numCols()-1)/2;
		SimpleMatrix It = new SimpleMatrix(Nx, Ny);
		
				
		
		return null;
	}

	private static SimpleMatrix lpf(SimpleMatrix sigma_n, double d) {
		// TODO Auto-generated method stub
		return null;
	}

	private static SimpleMatrix approxI1_I0(SimpleMatrix m) {
		// TODO Auto-generated method stub
		return null;
	}

	public static SimpleMatrix[] em_ml_rice2D(SimpleMatrix In, int N, double[] Ws) {
		double prod = Ws[0]*Ws[1];
		SimpleMatrix Mask = makeOnes((int)Ws[0], (int)Ws[1]);
		for (int i=0; i<Mask.numRows(); i++){
			for (int j=0; j<Mask.numCols(); j++){
				Mask.set(i, j, 1/prod);
			}
		}
		
//		a_k=sqrt(sqrt(max(2.*filter2B(Mask,In.^2).^2-filter2B(Mask,In.^4),0)));
		SimpleMatrix a_k = filter2b(Mask,In.elementPower(In));
		a_k = multipleM(a_k, 2.0);
		a_k.elementPower(a_k).minus(filter2b(Mask,In.elementPower(4.0)));
		a_k = max(a_k,0.0);
		a_k.elementPower(0.5);
		
//		sigma_k2=0.5.*max(filter2B(Mask,In.^2)-a_k.^2,0.01);
		SimpleMatrix sigma_k2 = filter2b(Mask,In.elementPower(In));
		sigma_k2.minus(a_k.elementPower(2.0));
		sigma_k2 = max(sigma_k2,0.01);
		sigma_k2 = multipleM(sigma_k2, 0.5);
		
		for(int i=1; i<N; i++){
			a_k = filter2b(Mask, approxI1_I0(a_k.elementMult(In).elementDiv(sigma_k2))).elementMult(In);
			a_k = max(a_k,0);
			sigma_k2 = filter2b(Mask, absM(In).elementPower(2));
			sigma_k2 = multipleM(sigma_k2, 0.5);
			sigma_k2.minus(divideM(a_k.elementPower(2),2.0));
			sigma_k2 = max(sigma_k2, 0.01);
		}
			
//		a_k=sqrt(sqrt(max(2.*filter2B(Mask,In.^2).^2-filter2B(Mask,In.^4),0)));
//		sigma_k2=0.5.*max(filter2B(Mask,In.^2)-a_k.^2,0.01);
//			
//		%EM ALGORITHM
//		for ii=1:N
//				a_k=max(filter2B(Mask,approxI1_I0(a_k.*In./sigma_k2).*In),0);
//				sigma_k2=max(0.5.*filter2B(Mask,abs(In).^2)-a_k.^2./2,0.01);
//				%sigma_k2=abs(0.5.*localmean3DB(abs(In).^2,Ws)-a_k.^2./2);
//		end
//
//		
//		 Signal=a_k;
//		Sigma_n=sqrt(sigma_k2);
		
		SimpleMatrix[] res = new SimpleMatrix[2];
		res[0] = a_k;
		res[1] = sigma_k2.elementPower(0.5);
		
		return res;
	}
	

	private static SimpleMatrix max(SimpleMatrix m, double d) {
		for (int i = 0; i < m.numRows(); i++) {
			for (int j = 0; j < m.numCols(); j++) {
				if (m.get(i,j) < d) {
					m.set(i, j, d);
				}
			}
		}
		return m;
	}

	private static SimpleMatrix divideM(SimpleMatrix m, double d) {
		for (int i=0; i<m.numRows(); i++){
			for (int j=0; j<m.numCols(); j++){
				double val = m.get(i, j)/d;
				m.set(i, j, val);
			}
		}
		return m;
	}

	private static SimpleMatrix absM(SimpleMatrix m) {
		for (int i=0; i<m.numRows(); i++){
			for (int j=0; j<m.numCols(); j++){
				m.set(i, j, Math.abs(m.get(i, j)));
			}
		}
		return m;
	}

	private static SimpleMatrix multipleM(SimpleMatrix m, double v){
		for (int i=0; i<m.numRows(); i++){
			for (int j=0; j<m.numCols(); j++){
				double val = m.get(i, j)*v;
				m.set(i, j, val);
			}
		}
		return m;
	}
	
}
