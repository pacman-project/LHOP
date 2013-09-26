package com.vicos.hop;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.BufferedOutputStream;
import java.io.BufferedInputStream;
import java.io.FileOutputStream;
import java.io.OutputStream;

public class HOPImage {

	public static final int HOP_FORMAT_BLOB = 0;
	public static final int HOP_FORMAT_PROTOBUF = 1;
	public static final int HOP_FORMAT_COMPRESS = 16; 


	private static String cfg1 = "type = struct; init_size = -150; scale_limit = 100;"
			+ "layer1_threshold = 0.1; power_correction = 0.6; scale_sigma = 1.0;"
			+ "scale_factor = 0.707107";

	private static String cfg2 = "layer_index = 2; layer_contraction = 2.0;"
			+ "r_response_threshold = 0.3; g_response_threshold = 0.3; candidate_r_threshold_vf = 0.8;"
			+ "candidate_g_threshold_vf = 0.8; realization_ratio_threshold = 0.9";

	private static String cfg3 = "layer_index = 3; layer_contraction = 2.0;"
			+ "r_response_threshold = 0.3; g_response_threshold = 0.3; candidate_r_threshold_vf = 0.8;"
			+ "candidate_g_threshold_vf = 0.8; realization_ratio_threshold = 0.9";

	private static void loadLibrary(final String name) throws IOException {
		String version = HOPImage.class.getPackage().getImplementationVersion();

		if (version == null) {
			version = "";
		}

		version = version.replace('.', '_');

		File f = new File(System.getProperty("java.io.tmpdir"), "lib" + name
				+ ".so");
		boolean exists = f.isFile(); // check if it already exists

		int bits = 32;

		// check for 64-bit systems and use to appropriate library
		if (System.getProperty("os.arch").indexOf("64") != -1)
			bits = 64;

		InputStream in = new BufferedInputStream(HOPImage.class
				.getResourceAsStream("lib" + name + ".so"));

		try {
			OutputStream fout = new BufferedOutputStream(
					new FileOutputStream(f));
			byte[] bytes = new byte[1024];

			for (int n = 0; n != -1; n = in.read(bytes)) {
				fout.write(bytes, 0, n);
			}

			fout.close();
		} catch (IOException ioe) {
			if (!exists) {
				throw ioe;
			}
		}

		f.deleteOnExit();

		System.load(f.getAbsolutePath());
	}

	static {
		try {
			loadLibrary("jhop");
		} catch (IOException e) {
			throw new RuntimeException(e);
		}

	}

	private hop_image image;

	private hop_library library;

	private HOPLayer[][] layers = null;

	private String name;

	public HOPImage(String name, File imageFile, File libraryFile)
			throws IOException {

		this.name = name;

		image = JHOP.hop_read_image(imageFile.toString());

		JHOP.hop_save_image(image, "test.png");
		
		library = JHOP.hop_read_library(libraryFile.toString());

	}

	public int getImageWidth() {
		return image.get_width();
	}

	public int getImageHeight() {
		return image.get_height();
	}

	public int getScalesCount() {
		return layers.length;
	}

	public int getLayersCount() {
		return layers[0].length;
	}

	public HOPLibrary getPartsLibrary() {
		return null;
	}

	public String getName() {
		return name;
	}

	public String getChannel() {
		return "hop.image." + name;
	}

	public void process() throws Exception {

		long start = System.currentTimeMillis();

		hop_result[] results = JHOP.hop_inference(image, cfg1);

		System.out.println(String.format("Processing 1st layer: %d ms (%d scales)", System.currentTimeMillis() - start, results.length));

		start = System.currentTimeMillis();
		
		for (int i = 0; i < results.length; i++)
			JHOP.hop_inference(results[i], library, cfg2);

		System.out.println(String.format("Processing 2nd layer: %d ms (%d scales)", System.currentTimeMillis() - start, results.length));
		
		start = System.currentTimeMillis();

		for (int i = 0; i < results.length; i++)
			JHOP.hop_inference(results[i], library, cfg3);

		System.out.println(String.format("Processing 3rd layer: %d ms (%d scales)", System.currentTimeMillis() - start, results.length));

		for (int i = 0; i < results.length; i++)
			results[i].to_file(String.format("scale%d.dat", i), HOP_FORMAT_PROTOBUF | HOP_FORMAT_COMPRESS);

	}

	public HOPLayer get(int scale, int layer) {
		return layers[scale][layer];
	}

	public static void main(String[] args) throws Exception {

		HOPImage image = new HOPImage("test", new File(
				"/home/lukacu/workspace/libhop/res/sample1.png"), new File(
				"/home/lukacu/workspace/libhop/res/layer5lib.plb"));

		image.process();
	}

}
