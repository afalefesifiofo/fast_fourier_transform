import java.io.IOException;


public class FFT_2D {

	//renvoie la TFD d'une image de complexes
	public static CpxImg FFT(CpxImg I) {

		CpxImg out = new CpxImg(I.taille());

		// FFT 1D sur les lignes
		for (int k = 0; k < I.taille(); k++)
			out.set_line(k,FFT_1D.FFT(I.get_line(k)));
		  
		// transposition
		out.transpose();

		// FFT 1D sur les "nouvelles" lignes de out (anciennes colonnes)
		for (int k = 0; k < I.taille(); k++)
			out.set_line(k,FFT_1D.FFT(out.get_line(k)));

		//on re transpose pour revenir dans le sens de d�part
		out.transpose();
		
		//on divise par la taille de I
		out.multiply(1./I.taille());
		return out.recentrage();
	}
	
	//renvoie la TFD inverse d'une images de complexes
	public static CpxImg FFT_inverse(CpxImg I) {
		I = I.recentrage();
		CpxImg out = new CpxImg(I.taille());
		for (int k = 0; k < I.taille(); k++)
			out.set_line(k, I.get_line(k).conjugue());

		out = FFT(out).recentrage();
		for (int k = 0; k < I.taille(); k++)
			out.set_line(k, out.get_line(k).conjugue());
		return out;
	}

	// compression par mise � z�ro des coefficients de fr�quence 
	// FI contient la TDF de I 
	// Dans FI on met � z�ros tous les coefficients correspondant � des fr�quences inf�rieures � k
	public static void compression(CpxImg FI, int k) {
		// A COMPLETER
	}

	// compression par seuillage des coefficients faibles
	// FI contient la TDF de I 
	// Dans FI on met � z�ros tous les coefficients dont le module est inf�rieur � seuil 
	// on renvoie le nombre de coefficients conserv�s 
	public static int compression_seuil(CpxImg FI, double seuil){
		//A COMPLETER
		return 0;
	}

	
	public static void main(String[] args) {
		
		try {			
			//PLACEZ ICI VOS TESTS en 2D
			//Exemple, lecture
			BytePixmap BP = new BytePixmap("AC_tp2_part2_donnees/tigre_512.pgm");
			CpxImg I = new CpxImg(BP);

			CpxImg tigrefft = FFT(I);
			BytePixmap tigreBP = tigrefft.convert_to_BytePixmap();

			tigreBP.write("tigrefft.pgm");

			CpxImg tigre_after_FFT_inverse = FFT_inverse(tigrefft);
			BytePixmap tigre_after_FFT_inverse_bp = tigre_after_FFT_inverse.convert_to_BytePixmap();

			tigre_after_FFT_inverse_bp.write("tigre_after_fft_inverse.pgm");

			// image compression
			System.out.println("Image compression\n");

			int n = tigrefft.taille();

			tigrefft.set_p_reel(n/2, n/2, 0);
			tigrefft.set_p_imag(n/2, n/2, 0); 

			CpxImg tigre_after_FFT_inverse2 = FFT_inverse(tigrefft);
			BytePixmap tigre_after_FFT_inverse_bp2 = tigre_after_FFT_inverse2.convert_to_BytePixmap();

			tigre_after_FFT_inverse_bp2.write("tigre_after_fft_inverse_compressed.pgm");

		} catch (IOException e) {
			e.printStackTrace();
		}
	}

}
