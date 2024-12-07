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
			int n = FI.taille(); // Taille du tableau (n x n)
		// Coordonnées du centre
		int centre = n / 2;
		// Parcourir tous les coefficients du tableau
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				// Vérifier si le coefficient est en dehors du carré centré de côté 2k
				if (Math.abs(i - centre) > k || Math.abs(j - centre) > k) {
					// Annuler le coefficient (mettre à 0)
					FI.set_p_reel(i, j, 0);
					FI.set_p_imag(i, j, 0);
				}
			}
		}
	}

	// compression par seuillage des coefficients faibles
	// FI contient la TDF de I 
	// Dans FI on met � z�ros tous les coefficients dont le module est inf�rieur � seuil 
	// on renvoie le nombre de coefficients conserv�s 
	public static int compression_seuil(CpxImg FI, double seuil){
		//A COMPLETER
		int n = FI.taille(); // Taille de l'image (n x n)
		int nbSignificatifs = 0; // Compteur de coefficients significatifs
	
		// Parcourir chaque coefficient de la matrice FI
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				// Obtenir les parties réelle et imaginaire
				double reel = FI.get_p_reel(i, j);
				double imaginaire = FI.get_p_imag(i, j);
	
				// Calculer le module du coefficient
				double module = Math.sqrt(reel * reel + imaginaire * imaginaire);
	
				// Vérifier si le module est significatif
				if (module >= seuil) {
					nbSignificatifs++; // Le coefficient est conservé
				} else {
					// Annuler les coefficients insignifiants
					FI.set_p_reel(i, j, 0);
					FI.set_p_imag(i, j, 0);
				}
			}
		}
		// Retourner le nombre de coefficients significatifs
		return nbSignificatifs;
	}

	
	public static void main(String[] args) {
		
		try {			
			//PLACEZ ICI VOS TESTS en 2D
			//Exemple, lecture
			BytePixmap BP = new BytePixmap("AC_tp2_part2_donnees/tigre_512.pgm");
			CpxImg I = new CpxImg(BP);
			CpxImg tigrefft = FFT(I);
			compression(tigrefft, 100);
			CpxImg tigre_after_FFT_inverseC = FFT_inverse(tigrefft);
			//BytePixmap tigreBP = tigrefft.convert_to_BytePixmap();
			//tigreBP.write("tigrefft.pgm");
			//CpxImg tigre_after_FFT_inverse = FFT_inverse(tigrefft);
			//BytePixmap tigre_after_FFT_inverse_bp = tigre_after_FFT_inverse.convert_to_BytePixmap();
			BytePixmap tigre_after_compression_bp = tigre_after_FFT_inverseC.convert_to_BytePixmap();

			tigre_after_compression_bp.write("tigre_after_compressionC.pgm");

			// image compression
			System.out.println("Image compression\n");

			int n = tigrefft.taille();

			tigrefft.set_p_reel(n/2, n/2, 0);
			tigrefft.set_p_imag(n/2, n/2, 0); 

			//CpxImg tigre_after_FFT_inverse2 = FFT_inverse(tigrefft);
			//BytePixmap tigre_after_FFT_inverse_bp2 = tigre_after_FFT_inverse2.convert_to_BytePixmap();

			//tigre_after_FFT_inverse_bp2.write("tigre_after_fft_inverse_compressed.pgm");

		} catch (IOException e) {
			e.printStackTrace();
		}
	}

}
