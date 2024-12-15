import java.io.IOException;


public class FFT_2D {

	public static CpxImg annulerCoefficient(CpxImg FI, int x, int y) {
		FI.set_p_reel(x, y, 0);
		FI.set_p_imag(x, y, 0);
		return FI;
	}	

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
		int centre = n / 2;
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				if (Math.abs(i - centre) > k || Math.abs(j - centre) > k) {
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
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				double reel = FI.get_p_reel(i, j);
				double imaginaire = FI.get_p_imag(i, j);
				double module = Math.sqrt(reel * reel + imaginaire * imaginaire);
				if (module >= seuil) {
					nbSignificatifs++; // Le coefficient est conservé
				} else {
					FI.set_p_reel(i, j, 0);
					FI.set_p_imag(i, j, 0);
				}
			}
		}
		return nbSignificatifs;
	}

	
	public static void main(String[] args) {
		
		try {			
			//PLACEZ ICI VOS TESTS en 2D
			//Exemple, lecture
			BytePixmap BP = new BytePixmap("AC_tp2_part2_donnees/barbara_512.pgm");
			CpxImg I = new CpxImg(BP);
			CpxImg tigrefft = FFT(I);
			String cheminImage = "AC_tp2_part2_donnees/barbara_512.pgm"; // Chemin de l'image
			String sortieCentral = "Resultats/tigre_sans_DC.pgm";
			String sortieBasseFreq = "Resultats/tigre_sans_basse.pgm";
			String sortieHauteFreq = "Resultats/tigre_sans_haute.pgm";
			int n = tigrefft.taille();
			int centre = n / 2;
			

			/*// Cas 1 : Annuler le coefficient central
			CpxImg result = annulerCoefficient(tigrefft, centre, centre);
			result = FFT_inverse(tfdImage);
			BytePixmap resultatSansDC = result.convert_to_BytePixmap();
			resultatSansDC.write(sortieCentral);
			System.out.println("Image sans DC sauvegardée sous : " + sortieCentral);

			// Cas 2 : Annuler un coefficient basse fréquence
			CpxImg result = annulerCoefficient(tigrefft,(tigrefft.taille()/2), tigrefft.taille()/2);
			result = FFT_inverse(result);
			BytePixmap resultatSansBF = result.convert_to_BytePixmap();
			resultatSansBF.write(sortieBasseFreq);
			System.out.println("Image sans basse fréquence sauvegardée sous : " + sortieBasseFreq);

		/*	// Cas 3 : Annuler un coefficient haute fréquence
			CpxImg tfdSansHauteFreq = tfdImage; // Copie
			CpxImg result = annulerCoefficient(tfdSansHauteFreq, n-1, n-1);
			BytePixmap resultatSansHauteFreq = FFT_inverse(result).convert_to_BytePixmap();
			resultatSansHauteFreq.write(sortieHauteFreq);
			System.out.println("Image sans haute fréquence sauvegardée sous : " + sortieHauteFreq);
*/              int nbCoeffs = (int) Math.round(0.1 * n * n); // Calcul dynamique
				System.out.println(n);
				compression_seuil(tigrefft, 40); // Appel de la fonction avec nbCoeffs
				//tigrefft.convert_to_BytePixmap().write("Resultats/fft_tigre_512.pgm");
				FFT_inverse(tigrefft).convert_to_BytePixmap().write("Resultats/compression_seuil_barbara.pgm");;
				BytePixmap tigreBP = tigrefft.convert_to_BytePixmap();
				 //image compression
				System.out.println("Image compression\n");
/* 
					n = tigrefft.taille();

				tigrefft.set_p_reel(n/2, n/2, 0);
				tigrefft.set_p_imag(n/2, n/2, 0); 
	*/
			//CpxImg tigre_after_FFT_inverse2 = FFT_inverse(tigrefft);
				//BytePixmap tigre_after_FFT_inverse_bp2 = tigre_after_FFT_inverse2.convert_to_BytePixmap();

				//tigre_after_FFT_inverse_bp2.write("tigre_after_fft_inverse_compressed.pgm");

		} catch (IOException e) {
			e.printStackTrace();
		}
	}

}
