package fr.cea.ig.metatarget.utils;

import java.io.File;
import java.io.IOException;
import java.net.URISyntaxException;
import java.net.URL;
//import java.nio.file.Files;
import java.text.DateFormat;
import java.text.NumberFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

//import org.apache.commons.lang3.RandomStringUtils;

//import org.mapdb.DB;
//import org.mapdb.DBMaker;

public class Utils {
	
	private static final DateFormat dateFormat = new SimpleDateFormat("yyyy/MM/dd HH:mm:ss");
	
	/*
	private static DB db = null;
	
	static{
		db = DBMaker
				.newHeapDB()
			    .transactionDisable()
			    .asyncWriteFlushDelay(100)
			    .make();
		//db = DBMaker.newMemoryDB().asyncWriteEnable(). /*cacheHardRefEnable()./transactionDisable().make();
		//db = DBMaker.newAppendFileDB(new File("/env/cns/bigtmp1/agkanogi/mapdb.temp")).closeOnJvmShutdown().transactionDisable().make();
		//db = DBMaker.newFileDB(new File("/env/cns/bigtmp1/agkanogi/mapdb.temp")).closeOnJvmShutdown().randomAccessFileEnable().transactionDisable().make();
	}
	
	public static DB getDB(){
		return db;
	}
	*/
		
	public static String time() {
		Calendar cal = Calendar.getInstance();
		return dateFormat.format(cal.getTime());
	}
	
	public static String RAMInfo(Runtime runtime){
		double gb = (double)1024*1024*1024;
        NumberFormat format = NumberFormat.getInstance(); 
        StringBuilder sb = new StringBuilder();
        
        sb.append("\n##### Heap utilization statistics [GB] #####\n");
        sb.append("Max Memory=" + format.format(runtime.maxMemory() / gb) + "\n");
        sb.append("Current Total Memory=" + format.format(runtime.totalMemory() / gb) + "\n");
        sb.append("Current Used Memory="  + format.format((runtime.totalMemory() - runtime.freeMemory()) / gb) + "\n");
        sb.append("Current Free Memory="  + format.format(runtime.freeMemory() / gb) + "\n");
        sb.append("############################################\n");
	
        return sb.toString();
	}

	final public static byte[][] concatSeqQual(final ArrayList<byte[]> arrays) {
		byte[] all = concat(arrays);
		int size = all.length;
		int sizeSeq = 0;
		int sizeQual = 0;
		for(byte b : all){
			if('+' == (char)(b & 0xff)){
				sizeQual = size - sizeSeq - 1;
				break;
			}
			else{
				sizeSeq++;
			}
		}
		
		byte[][] ret = new byte[2][];
		if(sizeSeq>0){
			byte[] destSeq = new byte[sizeSeq];
			System.arraycopy(all, 0, destSeq, 0, sizeSeq);
			ret[0] = destSeq;
		}
		if(sizeQual>0){
			byte[] destQual = new byte[sizeQual];
			System.arraycopy(all, sizeSeq+1, destQual, 0, sizeQual);
			ret[1] = destQual;
		}
		
		return ret;
	}
	
	final public static byte[] concat(final ArrayList<byte[]> arrays) {
		int size = 0;
		for (byte[] a : arrays)
			size += a.length;

		byte[] dest = new byte[size];

		int destPos = 0;
		for (int i = 0; i < arrays.size(); i++) {
			if (i > 0)
				destPos += arrays.get(i - 1).length;
			int length = arrays.get(i).length;
			System.arraycopy(arrays.get(i), 0, dest, destPos, length);
		}

		return dest;
	}
	
	public static int round(double d){
	    double dAbs = Math.abs(d);
	    int i = (int) dAbs;
	    double result = dAbs - (double) i;
	    if(result<0.5){
	        return d<0 ? -i : i;            
	    }else{
	        return d<0 ? -(i+1) : i+1;          
	    }
	}

	public static char[] bytesToStringUTFCustom(byte[] bytes) {
		char[] buffer = new char[bytes.length >> 1];
		for (int i = 0; i < buffer.length; i++) {
			int bpos = i << 1;
			char c = (char) (((bytes[bpos] & 0x00FF) << 8) + (bytes[bpos + 1] & 0x00FF));
			buffer[i] = c;
		}
		return buffer;
		//return new String(buffer);
	}
	
	public static byte[] stringToBytesUTFCustom(char[] str) {
		byte[] b = new byte[str.length << 1];
		for (int i = 0; i < str.length; i++) {
			char strChar = str[i];
			int bpos = i << 1;
			b[bpos] = (byte) ((strChar & 0xFF00) >> 8);
			b[bpos + 1] = (byte) (strChar & 0x00FF);
		}
		return b;
	}
	
	public static List<String> directoriesToFiles(List<String> inputFastaFileNames){
		ArrayList<String> outputFastaFileNames = new ArrayList<String>();
		try {
			for(String s : inputFastaFileNames){
				File f = new File(s);
				if(!f.isDirectory()){
					outputFastaFileNames.add(f.getCanonicalPath());
				}
				else{
					for(File ff : f.listFiles()){
						outputFastaFileNames.add(ff.getCanonicalPath());
					}
				}
			}
		} 
		catch (IOException e) {
			e.printStackTrace();
		}
		return outputFastaFileNames;
	}
	
	public static <K, V extends Comparable<? super V>> Map<K, V> sortByValue( Map<K, V> map ){
		List<Map.Entry<K, V>> list = new LinkedList<Map.Entry<K, V>>( map.entrySet() );
		Collections.sort( list, new Comparator<Map.Entry<K, V>>(){
			public int compare( Map.Entry<K, V> o1, Map.Entry<K, V> o2 ){
				int compare = (o1.getValue()).compareTo(o2.getValue());
				if(compare == 0){
					Long key1 = (Long)o1.getKey();
					Long key2 = (Long)o2.getKey();
					return +1 * (key1).compareTo(key2);
				}
				else{
					return -1 * compare;
				}
			}
		} );

		Map<K, V> result = new LinkedHashMap<K, V>();
		for (Map.Entry<K, V> entry : list){
			result.put( entry.getKey(), entry.getValue() );
		}
		return result;
	}
	
	public static ArrayList<String> generateKmerVocabulary(int k){
		ArrayList<String> vocabulary = new ArrayList<String>();
		
		char[] chars = "ATCG".toCharArray();
        int len = k;
        iterate(chars, len, new char[len], 0, vocabulary);
        
		return vocabulary;
	}
	
	private static void iterate(char[] chars, int len, char[] build, int pos, ArrayList<String> vocabulary) {
        if (pos == len) {
            String word = new String(build);
            vocabulary.add(word);
            return;
        }

        for (int i = 0; i < chars.length; i++) {
            build[pos] = chars[i];
            iterate(chars, len, build, pos + 1, vocabulary);
        }
    }
	
	/*
	public static String getRandomSequence(int length){ return
	RandomStringUtils.random(length, "ATCG"); }
	*/
	
	public static String reverseComplement(String sequence){
		StringBuilder sb = new StringBuilder();
		for(int i=sequence.length()-1; i>=0; i--){
			switch (sequence.charAt(i)){
				case 'A':
				case 'a':
					sb.append("T");
				break;

				case 'T':
				case 't':
					sb.append("A");
				break;

				case 'C':
				case 'c':
					sb.append("G");
				break;

				case 'G':
				case 'g':
					sb.append("C");
				break;

				default:
					;
				break;
			}
		}
		return sb.toString();
	}

	@SuppressWarnings("rawtypes")
	public static long classBuildTimeMillis(Class c){
	    URL resource = c.getResource(c.getSimpleName() + ".class");
	    if (resource == null) {
	        System.out.println("Failed to find class file for class: " + c.getName());
	        return 0L;
	    }
	    
	    if (resource.getProtocol().equals("file")) {
	    	try {
				return new File(resource.toURI()).lastModified();
			} 
	    	catch (URISyntaxException e) {
				e.printStackTrace();
				return 0L;
	    	}

	    } 
	    else if (resource.getProtocol().equals("jar")) {
	        String path = resource.getPath();
	        return new File(path.substring(5, path.indexOf("!"))).lastModified();
	    } 
	    else {
	        System.out.println("Unhandled url protocol: " + resource.getProtocol() + " for class: " + c.getName() + " resource: " + resource.toString());
	        return 0L;
	    }
	}

	public static <T extends Comparable<? super T>> List<T> asSortedList(Collection<T> c) {
		List<T> list = new ArrayList<T>(c);
		Collections.sort(list);
		return list;
	}
	
}
