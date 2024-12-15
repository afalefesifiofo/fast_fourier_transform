import java.awt.*;
import javax.swing.*;
import java.awt.geom.*;
import java.util.List;

public class Graph extends JPanel {
    private static final long serialVersionUID = 1L;
    List<Integer> tailles;
    List<Double> tempsFFT;
    List<Double> tempsCoef;

    int maxX;
    int maxY;
    JFrame frame;
    static final int marginX = 100;
    static final int marginY = 50;
    static final int arrowSize = 7;
    static final int radius = 5;
    static final int defaultWidth = 600;
    static final int defaultHeight = 600;
    static final int curveThickness = 2;

    public Graph(List<Integer> _tailles, List<Double> _tempsFFT, List<Double> _tempsCoef) {
        assert (_tailles.size() == _tempsFFT.size() && _tempsFFT.size() == _tempsCoef.size()) :
                "Les tableaux n'ont pas la même taille.";

        tailles = _tailles;
        tempsFFT = _tempsFFT;
        tempsCoef = _tempsCoef;
        // Trouver les maximums pour ajuster l'échelle
        maxX = tailles.stream().max(Integer::compareTo).orElse(1);
        maxY = Math.max(
                (int) Math.ceil(tempsFFT.stream().max(Double::compareTo).orElse(1.0)),
                (int) Math.ceil(tempsCoef.stream().max(Double::compareTo).orElse(1.0))
        );

        frame = new JFrame();
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
    }

    protected void paintComponent(Graphics g) {
        super.paintComponent(g);
        Graphics2D g2d = (Graphics2D) g;

        g2d.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
        int width = getWidth();
        int height = getHeight();

        // Dessiner les axes
        drawAxis(g2d, marginX, height - marginY, width - marginX, height - marginY, maxX, "Taille");
        drawAxis(g2d, marginX, height - marginY, marginX, marginY, maxY, "Temps");

        double scaleX = (double) (width - 2 * marginX) / maxX;
        double scaleY = (double) (height - 2 * marginY) / maxY;

        // Tracer la courbe FFT
        g2d.setPaint(Color.BLUE);
        g2d.setStroke(new BasicStroke(curveThickness));
        for (int i = 1; i < tailles.size(); i++) {
            double x1 = marginX + scaleX * tailles.get(i - 1);
            double y1 = height - marginY - scaleY * tempsFFT.get(i - 1);
            double x2 = marginX + scaleX * tailles.get(i);
            double y2 = height - marginY - scaleY * tempsFFT.get(i);
            g2d.draw(new Line2D.Double(x1, y1, x2, y2));
        }

        // Tracer la courbe Coefficients
        g2d.setPaint(Color.RED);
        g2d.setStroke(new BasicStroke(curveThickness));
        for (int i = 1; i < tailles.size(); i++) {
            double x1 = marginX + scaleX * tailles.get(i - 1);
            double y1 = height - marginY - scaleY * tempsCoef.get(i - 1);
            double x2 = marginX + scaleX * tailles.get(i);
            double y2 = height - marginY - scaleY * tempsCoef.get(i);
            g2d.draw(new Line2D.Double(x1, y1, x2, y2));
        }

        // Ajouter une légende
        g2d.setFont(new Font("Arial", Font.BOLD, 14));
        g2d.setColor(Color.BLUE);
        g2d.drawString("FFT", width - marginX - 100, marginY + 20);
        g2d.setColor(Color.RED);
        g2d.drawString("Coefficients", width - marginX - 100, marginY + 40);
    }

    private static void drawAxis(Graphics2D g, int x1, int y1, int x2, int y2, int maxValue, String legend) {
		Line2D.Double axis = new Line2D.Double(x1, y1, x2, y2);
		g.draw(axis);
	
		int numTicks = 10; // Ajustez pour avoir plus de ticks si nécessaire
		double step = (double) maxValue / numTicks;
	
		for (int i = 0; i <= numTicks; i++) {
			double pos = step * i;
			if (x1 == x2) { // Axe Y
				int yTick = (int) (y1 - (y1 - y2) * pos / maxValue);
				g.draw(new Line2D.Double(x1 - 5, yTick, x1 + 5, yTick));
				g.drawString(String.format("%.1f", pos), x1 - 35, yTick + 5);
			} else { // Axe X
				int xTick = (int) (x1 + (x2 - x1) * pos / maxValue);
				g.draw(new Line2D.Double(xTick, y1 - 5, xTick, y1 + 5));
				g.drawString(String.format("%.1f", pos), xTick - 10, y1 + 30);
			}
		}
	
		g.setFont(new Font("Arial", Font.PLAIN, 16));
		if (x1 == x2) { // Axe Y
			g.drawString(legend, x1 - 100, (y1 + y2) / 2); // Ajustez ici (précédemment -70)
		} else { // Axe X
			g.drawString(legend, (x1 + x2) / 2, y1 + 60); // Ajustez ici (précédemment +40)
		}
	}
	
    public void display() {
		frame.add(this);
		frame.setSize(800 + 2 * marginX, defaultHeight + 2 * marginY); // Largeur augmentée
		frame.setLocation(200, 200);
		frame.setVisible(true);
	}
}