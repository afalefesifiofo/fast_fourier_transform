/********************************************************
/*   La classe CpxTab représente une image carré 
 * (n x n) dont les éléments sont des nombres complexes
 ********************************************************/

public class CpxImg {


	/******** ATTRIBUT **********/
	// lignes[k] est la kème ligne de l'image
	private CpxTab[] lignes;
	
	/******** CONSTRUCTEURS **********/
	//construit une image de taille nxn (de valeurs nulles)
	public CpxImg(int n) {
		lignes = new CpxTab[n];
		for(int k=0;k<n;k++)
			lignes[k] = new CpxTab(n);
	}
	
	//construit un CpxImg à partir d'un BytePixmap
	public CpxImg(BytePixmap p) {
		this(p.height);
		assert (p.height==p.width) : "CpxImg à partir de BytePixmap : p.height!=p.width";

		double[] t = p.getDoubles();
		for (int k = 0; k < this.taille(); k++) 
			for (int l = 0; l < this.taille(); l++){
				this.set_p_reel(k, l, t[k * this.taille() + l]);
			}
	}
	
	/******** ACCESSEURS **********/
	//renvoie une dimension de l'image
	int taille(){
		return lignes.length;
	}
	//renvoie la partie réelle du pixel d'indice (k,l) 
	double get_p_reel(int k, int l){
		return lignes[k].get_p_reel(l);
	}
	//renvoie la partie imaginaire du pixel d'indice (k,l)
	double get_p_imag(int k, int l){
		return lignes[k].get_p_imag(l);
	}
	
	//renvoie la k-ième ligne
	CpxTab get_line(int k){
		return lignes[k];
	}
	
	//modifie la partie réelle du  k-ième complexes
	void set_p_reel(int k, int l, double x){
		lignes[k].set_p_reel(l, x);
	}
	//modifie la partie imaginaire du  k-ième complexes
	void set_p_imag(int k, int l, double x){
		lignes[k].set_p_imag(l, x);
	}
	
	//modifie la k-ième ligne
	void set_line(int k, CpxTab line){
		lignes[k] = line;
	}
	
	
	/******** METHODES **********/
	//transpose l'image : I(k,l) <- I(l,k)
	public void transpose(){
		for (int k = 0; k < this.taille(); k++) {
			for (int l = k + 1; l < this.taille(); l++) {
				double temp_r = this.get_p_reel(k, l);
				double temp_i = this.get_p_imag(k, l);
				this.set_p_reel(k, l, this.get_p_reel(l, k));
				this.set_p_imag(k, l, this.get_p_imag(l, k));
				this.set_p_reel(l, k, temp_r);
				this.set_p_imag(l, k, temp_i);
			}
		}
	}
	
	//multiplie chaque élément par x
	public void multiply(double x){
		for (int k = 0; k < this.taille(); k++) {
			for (int l = 0; l < this.taille(); l++) {
				this.set_p_reel(k, l, x*this.get_p_reel(k, l));
				this.set_p_imag(k, l, x*this.get_p_imag(k, l));
			}
		}
	}
	
	//recentre l'image pour l'afficher dans l'ordre classique dans l'espace de Fourier 2D 
	public CpxImg recentrage() {
		CpxImg out = new CpxImg(this.taille());
		int n = this.taille();
		for (int i = 0; i < n / 2; i++) {
			for (int j = 0; j < n; j++) {
				int k = (i + n / 2) % n;
				int l = (j + n / 2) % n;
				out.set_p_reel(i, j, this.get_p_reel(k, l));
				out.set_p_imag(i, j, this.get_p_imag(k, l));
				out.set_p_reel(k, l, this.get_p_reel(i, j));
				out.set_p_imag(k, l, this.get_p_imag(i, j));
			}
		}
		return out;
	}
	
	
	//convertit l'image en un BytePixmap
	//en sauvant le module de chaque complexes
	public BytePixmap convert_to_BytePixmap() {
		int n = this.taille();
		byte[] b = new byte[n * n];
		for (int k = 0; k < n; k++) {
			for (int l = 0; l < n; l++) {
				double reel = this.get_p_reel(k, l);
				double img  = this.get_p_imag(k, l);
				double module = Math.sqrt(reel*reel + img*img);				
				byte module_b = ((byte) Math.max(0,Math.min(255,module)));
				b[k * n + l] = module_b;
				//b[k * n + l] = ((byte) Math.sqrt(reel * reel + img*img) > 255 ? (byte) 255 : (byte) Math.sqrt(reel*reel+img*img));
			}
		}
		return new BytePixmap(n, n, b);
	}	
}
