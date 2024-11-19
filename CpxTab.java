/**********************************************
/*   La classe CpxTab représente un tableau
 *              de nombres complexes
 **********************************************/

public class CpxTab {
	
	
		/******** ATTRIBUTS **********/
		//p_r[k] contient la partie réelle du k-ème complexe
		private double[] p_r;
		
		//p_i[k] contient la partie imaginaire du k-ème complexe
		private double[] p_i;

		
		/******** CONSTRUCTEURS **********/
		//Construit un tableau de complexes de taille n de valeur nulle
		public CpxTab(int n) {
			this.p_r = new double[n];
			this.p_i = new double[n];
		}

		
		//construit un tableau de complexes à partir d'un tableau de réels
		public CpxTab(double[] t) {
			this.p_r = t; // la partie réelle est directement le tableau t
			this.p_i = new double[t.length]; // la partie imaginaire est nulle
		}

		//clone le tableau de complexes t
		public CpxTab(CpxTab t) {
			this.p_r = t.p_r.clone();
			this.p_i = t.p_i.clone();
		}
		
		/******** ACCESSEURS **********/
		//renvoie la taille du vecteur
		int taille(){
			return p_r.length;
		}
		//renvoie la partie réelle du  k-ième complexes 
		double get_p_reel(int k){
			return p_r[k];
		}
		//renvoie la partie imaginaire du  k-ième complexes
		double get_p_imag(int k){
			return p_i[k];
		}
		//modifie la partie réelle du  k-ième complexes
		void set_p_reel(int k, double x){
			p_r[k] = x;
		}
		//modifie la partie imaginaire du  k-ième complexes
		void set_p_imag(int k,double x){
			p_i[k] = x;
		}
		
		/******** METHODES **********/
		//renvoie le tableau constitué des conjugués du tableau this
		// si this[k] = x+ i*y alors  out[k] = x - i*y 
		public CpxTab conjugue() {
			//on clone this
			CpxTab out = new CpxTab(this);
			for (int k = 0; k < taille(); k++) {
				out.p_i[k] = -out.p_i[k]; //on change le signe de la partie imaginaire
			}
			return out;
		}

		//multiplie les deux tableaux terme à terme
		// out[k] = c1[k]*c2[k] (multiplication complexe)
		// c1 et c2 doivent être de même taille
		public static CpxTab multiplie(CpxTab c1, CpxTab c2) {
			assert (c1.taille()==c2.taille()) :"multiplie: c1 et c2 ne sont pas de même taille, taille c1="+c1.taille()+" taille c2="+c2.taille();
			CpxTab out = new CpxTab(c1.taille());
			for (int k = 0; k < c1.taille(); k++) {
				out.p_r[k] = (c1.p_r[k] * c2.p_r[k] - c1.p_i[k] * c2.p_i[k]);
				out.p_i[k] = (c1.p_r[k] * c2.p_i[k] + c1.p_i[k] * c2.p_r[k]);
			}
			return out;
		}

		// méthodes d'affichage
		// pour avoir au plus 9 décimales
		private static final long PRECISION = 1000000000;

		private String print(double x) {
			return "" + ((double) Math.round(x * PRECISION)) / PRECISION;
		}

		private static final double THRESHOLD = 0.0001;

		public String toString() {
			String r = "[";
			for (int k = 0; k < taille(); k++) {
				if (Math.abs(p_r[k]) > THRESHOLD) {
					r = r + print(p_r[k]);
					if (p_i[k] > THRESHOLD) {
						r += "+" + print(p_i[k]) + "i";
					} else if (p_i[k] < -THRESHOLD) {
						r += "-" + print(-p_i[k]) + "i";
					}
				} else {
					
					if (Math.abs(p_i[k]) < THRESHOLD) {
						r = r + "0";
					} else {
						r = r + print(p_i[k]) + "i";
					}
				}
				if (k < taille() - 1)
					r = r + " ";
			}
			r = r + "]";
			return r;
		}
	}