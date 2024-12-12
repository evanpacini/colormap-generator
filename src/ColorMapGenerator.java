import java.awt.image.BufferedImage;

public class ColorMapGenerator {
    private static final double[][] palette = readPalette("rust-gold-8-1x.png");
    private static final int[] intPalette = toIntPalette();

    private static double[][] readPalette(String path) {
        try {
            BufferedImage img = javax.imageio.ImageIO.read(new java.io.File(path));
            double[][] palette = new double[img.getWidth()][3];
            for (int i = 0; i < img.getWidth(); i++) {
                int rgb = img.getRGB(i, 0);
                palette[i] = new double[]{((rgb >> 16) & 0xff) / 255., ((rgb >> 8) & 0xff) / 255., (rgb & 0xff) / 255.};
            }
            return palette;
        } catch (java.io.IOException e) {
            throw new RuntimeException("Error while reading palette", e);
        }
    }

    private static int[] toIntPalette() {
        int[] intPalette = new int[palette.length];
        for (int i = 0; i < palette.length; i++)
            intPalette[i] = ((int) Math.round(palette[i][0] * 255) << 16) | ((int) Math.round(palette[i][1] * 255) << 8) | (int) Math.round(palette[i][2] * 255);
        return intPalette;
    }

    private static double srgbToLinear(double c) {
        if (c <= 0.04045) return c / 12.92;
        return Math.pow((c + 0.055) / 1.055, 2.4);
    }

    private static double[] srgbToLinear(double[] c) {
        return new double[]{srgbToLinear(c[0]), srgbToLinear(c[1]), srgbToLinear(c[2])};
    }

    private static double[] linearToXyz(double[] c) {
        return new double[]{0.4124564 * c[0] + 0.3575761 * c[1] + 0.1804375 * c[2], 0.2126729 * c[0] + 0.7151522 * c[1] + 0.0721750 * c[2], 0.0193339 * c[0] + 0.1191920 * c[1] + 0.9503041 * c[2]};
    }

    private static double[] xyzToLab(double[] c) {
        double[] white = new double[]{0.95047, 1.00000, 1.08883};
        for (int i = 0; i < 3; i++) c[i] /= white[i];
        for (int i = 0; i < 3; i++) c[i] = c[i] > 0.008856 ? Math.pow(c[i], 1. / 3) : 7.787 * c[i] + 16. / 116;
        return new double[]{116 * c[1] - 16, 500 * (c[0] - c[1]), 200 * (c[1] - c[2])};
    }

    private static double CIEDE2000(double[] lab1, double[] lab2) {
        double C_7 = Math.pow((Math.sqrt(lab1[1] * lab1[1] + lab1[2] * lab1[2]) + Math.sqrt(lab2[1] * lab2[1] + lab2[2] * lab2[2])) / 2., 7.);
        double G = 1. + (1. - Math.sqrt(C_7 / (C_7 + 6103515625.))) / 2.;

        double a1p = G * lab1[1];
        double a2p = G * lab2[1];

        double C1p = Math.sqrt(a1p * a1p + lab1[2] * lab1[2]);
        double C2p = Math.sqrt(a2p * a2p + lab2[2] * lab2[2]);

        double CpProd = C1p * C2p;

        double h1p = Math.atan2(lab1[2], a1p);
        double h2p = Math.atan2(lab2[2], a2p);

        if (h1p < 0) h1p += 2 * Math.PI;
        if (h2p < 0) h2p += 2 * Math.PI;

        double dhp = 0.0;
        if (CpProd != 0) {
            dhp = h2p - h1p;
            if (dhp > Math.PI) dhp -= 2 * Math.PI;
            else if (dhp < -Math.PI) dhp += 2 * Math.PI;
        }

        double dH = 2 * Math.sqrt(CpProd) * Math.sin(dhp / 2);

        double Cp = (C1p + C2p) / 2;
        double hp = (h1p + h2p) / 2;

        if (Math.abs(h1p - h2p) > Math.PI) hp -= Math.PI;
        if (hp < 0) hp += 2 * Math.PI;

        double LpSqr = Math.pow((lab1[0] + lab2[0]) / 2 - 50, 2);

        double T = 1 - 0.17 * Math.cos(hp - 0.52359877559829887307710723054658) + 0.24 * Math.cos(2 * hp) + 0.32 * Math.cos(3 * hp + 0.10471975511965977461542144610932) - 0.20 * Math.cos(4 * hp - 1.0995574287564276334619251841478);

        double deltaThetaRad = 1.0471975511965977461542144610932 * Math.exp(-5.25249016001879 * Math.pow(hp - 4.799655442984406, 2.));

        double Cp7 = Math.pow(Cp, 7);

        double RT = -Math.sin(deltaThetaRad) * 2 * Math.sqrt(Cp7 / (Cp7 + 6103515625.0));

        double dCSc = (C2p - C1p) / (1 + 0.045 * Cp);
        double dHSh = dH / (1 + 0.015 * Cp * T);
        return Math.pow((lab2[0] - lab1[0]) / (1 + 0.015 * LpSqr / Math.sqrt(20 + LpSqr)), 2.) + Math.pow(dCSc, 2.) + Math.pow(dHSh, 2.) + RT * dCSc * dHSh;
    }

    public static void main(String[] args) {
        int[] pixels = new int[16777216];
        java.util.stream.IntStream.range(0, 256).parallel().forEach(r -> {
            for (int g = 0; g < 256; g++) {
                for (int b = 0; b < 256; b++) {
                    double[] color = new double[]{(double) r / 255., (double) g / 255., (double) b / 255.};
                    color = xyzToLab(linearToXyz(srgbToLinear(color)));
                    int nearest = 0;
                    double minDistance = CIEDE2000(xyzToLab(linearToXyz(srgbToLinear(palette[0]))), color);
                    for (int i = 1; i < palette.length; i++) {
                        double distance = CIEDE2000(xyzToLab(linearToXyz(srgbToLinear(palette[i]))), color);
                        if (distance < minDistance) {
                            minDistance = distance;
                            nearest = i;
                        }
                    }
                    pixels[(b << 16) | (g << 8) | r] = intPalette[nearest];
                }
            }
        });
        BufferedImage img = new BufferedImage(4096, 4096, BufferedImage.TYPE_INT_RGB);
        img.setRGB(0, 0, 4096, 4096, pixels, 0, 4096);
        try {
            javax.imageio.ImageIO.write(img, "png", new java.io.File("colormap.png"));
        } catch (java.io.IOException ignored) {
        }
    }
}

