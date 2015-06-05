package homomorphic_sense;

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
		SimpleMatrix lRn = new SimpleMatrix();
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

	private static SimpleMatrix correct_rice_gauss(SimpleMatrix sNR) {
		// TODO Auto-generated method stub
		return null;
	}

	private static SimpleMatrix filter2b(SimpleMatrix ones, SimpleMatrix in) {
		// TODO Auto-generated method stub
		return null;
	}

	private static SimpleMatrix lpf(SimpleMatrix sigma_n, double d) {
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
		a_k.elementPower(a_k).minus(filter2b(Mask,In.elementPower(4.0)));
		for (int i=0; i<a_k.numRows(); i++){
			for (int j=0; j<a_k.numCols(); j++){
				double val = a_k.get(i, j)*2;
				a_k.set(i, j, val);
			}
		}
		
		SimpleMatrix sigma_k2 = filter2b(Mask,In.elementPower(In));
		sigma_k2.elementPower(2).minus(a_k.elementPower(2.0));
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
		res[0] = new SimpleMatrix();
		res[1] = new SimpleMatrix();
		
		return res;
	}
	
}
