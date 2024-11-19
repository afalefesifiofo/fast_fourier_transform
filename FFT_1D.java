/**** on va ici implémenter la transformée de Fourier rapide 1D ****/

public class FFT_1D {
	
	//"combine" c1 et c2 selon la formule vue en TD
	// c1 et c2 sont de même taille
	// la taille du résultat est le double de la taille de c1
	public static CpxTab combine(CpxTab c1, CpxTab c2) {
		assert (c1.taille() == c2.taille()) : 
			"combine: c1 et c2 ne sont pas de même taille, taille c1=" + c1.taille() + " taille c2=" + c2.taille();
	
		int n = 2 * c1.taille();  // La taille du résultat est deux fois celle de c1 (et c2)
		CpxTab resultat = new CpxTab(n);  // Le tableau pour stocker le résultat
	
		double pi2 = 2 * Math.PI / n;  // Calcul de 2π/n
	
		for (int k = 0; k < c1.taille(); k++) {
			// Calcul de ω_n^k = cos(2πk/n) + i*sin(2πk/n)
			double cos = Math.cos(pi2 * k);
			double sin = Math.sin(pi2 * k);
	
			// Extraction des parties réelle et imaginaire de c1[k] et c2[k]
			double c1Re = c1.get_p_reel(k);
			double c1Im = c1.get_p_imag(k);
			double c2Re = c2.get_p_reel(k);
			double c2Im = c2.get_p_imag(k);
	
			// Calcul de ω_n^k * c2[k]
			double omegaC2Re = c2Re * cos - c2Im * sin;
			double omegaC2Im = c2Re * sin + c2Im * cos;
	
			// Calcul de A(ω_k^n) = c1[k] + ω_n^k * c2[k]
			double resRe1 = c1Re + omegaC2Re;
			double resIm1 = c1Im + omegaC2Im;
	
			// Calcul de A(ω_k+n/2^n) = c1[k] - ω_n^k * c2[k]
			double resRe2 = c1Re - omegaC2Re;
			double resIm2 = c1Im - omegaC2Im;
	
			// Enregistrement des résultats dans le tableau de sortie
			resultat.set_p_reel(k, resRe1);
			resultat.set_p_imag(k, resIm1);
			resultat.set_p_reel(k + c1.taille(), resRe2);
			resultat.set_p_imag(k + c1.taille(), resIm2);
		}
	
		return resultat;
	}
	
	
	
	
	
	
	
	

	//renvoie la TFD d'un tableau de complexes
	//la taille de x doit être une puissance de 2
	public static CpxTab FFT(CpxTab x) {
		int n = x.taille();  // Get the size of the input array x
	
		// Base case: If the size is 1, return the array as is
		if (n == 1) {
			return x;
		}
	
		// Check that the size is a power of 2
		assert (n % 2 == 0) : "FFT: The size of x must be a power of 2";
	
		// Create two new CpxTab arrays to store even and odd indexed elements
		CpxTab even = new CpxTab(n / 2);  // Array for even-indexed elements
		CpxTab odd = new CpxTab(n / 2);   // Array for odd-indexed elements
	
		// Populate the "even" and "odd" arrays
		for (int i = 0; i < n / 2; i++) {
			// Fill the "even" array with elements at even indices
			even.set_p_reel(i, x.get_p_reel(2 * i));  // Real part of even indexed elements
			even.set_p_imag(i, x.get_p_imag(2 * i));  // Imaginary part of even indexed elements
	
			// Fill the "odd" array with elements at odd indices
			odd.set_p_reel(i, x.get_p_reel(2 * i + 1)); // Real part of odd indexed elements
			odd.set_p_imag(i, x.get_p_imag(2 * i + 1)); // Imaginary part of odd indexed elements
		}
	
		// Recursively apply the FFT to the "even" and "odd" arrays
		CpxTab evenFFT = FFT(even);
		CpxTab oddFFT = FFT(odd);
	
		// Combine the results of the even and odd FFTs using the combine function
		return combine(evenFFT, oddFFT);
	}
	

	//renvoie la TFD d'un tableau de réels
	//la taille de x doit être une puissance de 2
	public static CpxTab FFT(double[] x) {
		return FFT(new CpxTab(x));
	}
			
	//renvoie la transformée de Fourier inverse de y
	public static CpxTab FFT_inverse(CpxTab y) {
		//A FAIRE
		return null;
	}
	
	//calcule le produit de deux polynômes en utilisant la FFT
	//tab1 et tab2, sont les coefficients de ces polynômes
	// CpxTab sera le tableau des coefficients du polynôme produit (purement réel)
	public static CpxTab multiplication_polynome_viaFFT(double[] tab1, double[] tab2) {
		
		//on commence par doubler la taille des vecteurs en rajoutant des zéros à la fin (cf TD)
		double[] t1 = new double[2*tab1.length], t2 = new double[2*tab1.length];
		for (int i = 0; i < tab1.length; i++) {
			t1[i] = tab1[i];
			t2[i] = tab2[i];
		}

		//A COMPLETER !!
		return null;
	}

	
	//renvoie un tableau de réels aléatoires
	//utile pour tester la multiplication de polynômes
	public static double[] random(int n) {
		double[] t = new double[n];

		for (int i = 0; i < n; i++)
			t[i] = Math.random();
		return t;
	}

	//effectue la multiplication de polynômes représentés par coefficients
	// p1, p2 les coefficients des deux polynômes P1 et P2
	// renvoie les coefficients du polynôme P1*P2
	private static double [] multiplication_polynome_viaCoeff(double[] p1, double[] p2){
		
		int n = p1.length + p2.length - 1;
		double a,b;
		double [] out = new double[n];
		for (int k = 0; k < n; k++) {
			for (int i = 0; i <= k; i++) {
				a = (i<p1.length) ? p1[i]:0;
				b = (k-i<p2.length) ? p2[k-i] : 0;
				out[k] += a*b;
			}
		}
		return out;
	}
	

	//affiche un tableau de réels
	private static void afficher(double [] t){
		System.out.print("[");
		for(int k=0;k<t.length;k++){
			System.out.print(t[k]);
			if (k<(t.length-1))
				System.out.print(" ");
		}
		System.out.println("]");
	}
	
	public static void main(String[] args) {
		double[] t5 = {1,2,3,4};

		 // Coefficienti del polinomio: [1, 2, 3, 4]
		 double[] realPart = {1, 2, 3, 4};
		 CpxTab input = new CpxTab(realPart);
 
		 // Dividiamo i coefficienti in pari e dispari
		 double[] even = {1, 3}; // coefficienti pari
		 double[] odd = {2, 4};  // coefficienti dispari
		 CpxTab c1 = new CpxTab(even); // TFD dei pari
		 CpxTab c2 = new CpxTab(odd);  // TFD dei dispari
 
		 // Calcolo della TFD combinata
		 CpxTab result = FFT_1D.combine(c1, c2);
 
		 // Stampiamo il risultato
		 System.out.println("Risultato TFD combinata: " + result);
 

		/* Exo 2: calculez et affichez TFD(1,2,3,4) */
		CpxTab x = new CpxTab(4);
		x.set_p_reel(0, 1.0);  // Real part of x[0]
		x.set_p_imag(0, 0.0);  // Imaginary part of x[0]
		x.set_p_reel(1, 2.0);  // Real part of x[1]
		x.set_p_imag(1, 0.0);  // Imaginary part of x[1]
		x.set_p_reel(2, 3.0);  // Real part of x[2]
		x.set_p_imag(2, 0.0);  // Imaginary part of x[2]
		x.set_p_reel(3, 4.0);  // Real part of x[3]
		x.set_p_imag(3, 0.0);  // Imaginary part of x[3]
	
		// Apply the FFT
		CpxTab result_fft = FFT_1D.FFT(x);

		System.out.println(result_fft);
		
		/* Exo 3: calculez et affichez TFD_inverse(TFD(1,2,3,4)) */
			//A FAIRE		

		/* Exo 4: multiplication polynomiale, vérification*/
			/* A(X) = 2 et B(X)=-3 */
			//A FAIRE		

			/* A(X) = 2+X et B(X)= -3+2X */
			//A FAIRE					

			/* A(X) = 1 + 2X + 3X^2 + 4X^3 et B(X) = -3 + 2X - 5 X^2*/
	/*
		System.out.println("-----------------------------------------------------");
		System.out.println("   Comparaison des 2 méthodes de multiplications polynomiales");
		double[] t6 = {-3,2,-5,0};
		System.out.println("mult via FFT  --> " + multiplication_polynome_viaFFT(t5, t6));
		System.out.print(  "mult via coeff -> ");
		afficher(multiplication_polynome_viaCoeff(t5, t6));
	*/

		/* Exo 5: comparaison des temps de calculs */
	/*
		// Pour étude du temps de calcul 
		int n = 256;  // taille des polynômes à multiplier (testez différentes valeurs en gardant des puissances de 2)
			
		System.out.println("Temps de calcul pour n="+n);
		double[] tab1 =random(n),tab2 = random(n);
		long date1, date2;
		date1 = System.currentTimeMillis();
		multiplication_polynome_viaCoeff(tab1, tab2);
		date2 = System.currentTimeMillis();
		System.out.println("   via Coeff: " + (date2 - date1));

		date1 = System.currentTimeMillis();
		multiplication_polynome_viaFFT(tab1, tab2);
		date2 = System.currentTimeMillis();
		System.out.println("   via FFT  : " + (date2 - date1));
	*/

	}

}