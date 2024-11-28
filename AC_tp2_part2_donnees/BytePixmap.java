import java.io.*;

public class BytePixmap extends Pixmap {

  public final byte[] data;

  public BytePixmap(int w, int h, byte[] pixels) throws IllegalArgumentException {
    super(w, h);
    if ( pixels == null )
      pixels = new byte[w*h];
    if ( pixels.length != w*h )
      throw new IllegalArgumentException("bad dimensions");
    data = pixels;
  }

  public BytePixmap(int w, int h) {
    this(w, h, null);
  }

  public BytePixmap(Pixmap p) {
    this(p.width, p.height, p.getBytes());
  }

  public BytePixmap(String fileName) throws IOException {
    super(fileName);
    data = readBytes();
    close();
  }

  // conversions
  public byte[] getBytes() {
    return getBytes(data);
  }

  public short[] getShorts() {
    return getShorts(data);
  }

  public double[] getDoubles() {
    return getDoubles(data);
  }

  // strings
  public String pixelType() {
    return "byte";
  }

}
