import java.io.*;

class PixmapReader extends FileInputStream {

  private char c;

  public PixmapReader(String fileName) throws FileNotFoundException {
    super(fileName);
  }

  public boolean matchKey(String key) throws IOException {
    byte[] buf = new byte[key.length()];
    read(buf, 0, key.length());
    return key.compareTo(new String(buf)) == 0;
  }

  public void getChar() throws IOException {
    c = (char)read();
  }

  public int getInt() throws IOException {
    String s = "";
    while ( (c != '\n') && Character.isSpaceChar(c) ) 
      getChar();
    while ( (c != '\n') && !Character.isSpaceChar(c) ) {
      s = s + c;
      getChar();
    }      
    return Integer.parseInt(s);
  }

  public void skipLine() throws IOException {
    while ( c != '\n' )
      getChar();
  }

  public void skipComment(char code) throws IOException {
    getChar();
    while ( c == code ) {
      skipLine();
      getChar();
    }
  }

  public byte[] loadData(int size) throws IOException {
    byte[] data = new byte[size];
    read(data, 0, size);
    return data;
  }

  public void close() {
    try {
      super.close();
    } catch (IOException e) {}
  }

}

class PixmapWriter extends FileOutputStream {

  public PixmapWriter(String fileName) throws IOException {
    super(fileName);
  }

  public void put(String data) throws IOException {
    write(data.getBytes());
  }

}

public abstract class Pixmap {

  private static String MAGIC_PGM = "P5\n";

  public final int width;
  public final int height;
  public final int size;
  private PixmapReader reader;

  public Pixmap(int w, int h) {
    width = w;
    height = h;
    size = width * height;
  }

  Pixmap(String fileName, String magic) throws IOException {
    reader = new PixmapReader(fileName);
    if ( !reader.matchKey(magic) )
      throw new IOException(fileName + " : wrong magic number");
    reader.skipComment('#');
    width = reader.getInt();
    height = reader.getInt();
    size = width * height;
    reader.skipLine();
    reader.skipComment('#');
    reader.skipLine();
  }

  public Pixmap(String fileName) throws IOException {
    this(fileName, MAGIC_PGM);
  }

  final byte[] readBytes(int size) throws IOException {
    return reader.loadData(size);
  }

  final byte[] readBytes() throws IOException {
    return readBytes(size);
  }

  final void close() {
    reader.close();
    reader = null;
  }

  final void write(String fileName, String magic, byte[] buffer) {
    try {
      PixmapWriter writer = new PixmapWriter(fileName);
      writer.put(magic);
      writer.put("#\n");
      writer.put(width + " " + height + "\n255\n");
      writer.write(buffer);
      writer.close();
      System.err.println("'" + fileName + "' wrote");
    } catch (IOException e) {
      System.err.println("can't write '"+fileName+"'");
    }
  }

  public void write(String fileName) {
    write(fileName, MAGIC_PGM, getBytes());
  }

  // conversions
  public abstract byte[] getBytes();
  public abstract short[] getShorts();
  public abstract double[] getDoubles();

  public static int intValue(byte b) {
    if ( b < 0 )
      return b + 256;
    else
      return b;
  }

  public static byte[] getBytes(byte[] buffer) {
    byte[] data = new byte[buffer.length];
    for (int i = 0; i < buffer.length; i++)
      data[i] = buffer[i];
    return data;
  }

  public static byte[] getBytes(short[] buffer) {
    byte[] data = new byte[buffer.length];
    for (int i = 0; i < buffer.length; i++)
      data[i] = (byte)buffer[i];
    return data;
  }

  public static byte[] getBytes(double[] buffer) {
    byte[] data = new byte[buffer.length];
    for (int i = 0; i < buffer.length; i++)
      data[i] = (byte)buffer[i];
    return data;
  }

  public static short[] getShorts(byte[] buffer) {
    short[] data = new short[buffer.length];
    for (int i = 0; i < buffer.length; i++)
      data[i] = (short)intValue(buffer[i]);
    return data;
  }

  public static short[] getShorts(short[] buffer) {
    short[] data = new short[buffer.length];
    for (int i = 0; i < buffer.length; i++)
      data[i] = buffer[i];
    return data;
  }

  public static short[] getShorts(double[] buffer) {
    short[] data = new short[buffer.length];
    for (int i = 0; i < buffer.length; i++)
      data[i] = (short)buffer[i];
    return data;
  }

  public static double[] getDoubles(byte[] buffer) {
    double[] data = new double[buffer.length];
    for (int i = 0; i < buffer.length; i++)
      data[i] = intValue(buffer[i]);
    return data;
  }

  public static double[] getDoubles(short[] buffer) {
    double[] data = new double[buffer.length];
    for (int i = 0; i < buffer.length; i++)
      data[i] = buffer[i];
    return data;
  }

  public static double[] getDoubles(double[] buffer) {
    double[] data = new double[buffer.length];
    for (int i = 0; i < buffer.length; i++)
      data[i] = buffer[i];
    return data;
  }

  // strings

  public abstract String pixelType();

  public String toString() {
    return "["+pixelType()+"-map : W = "+width+"  H = "+height+"]";
  }

}
