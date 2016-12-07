package by.vma.lab8;

public class Main {
    private static class Matrix {
        public double[][] matrix;
        private int lines;
        private int columns;

        public Matrix(int lines, int columns) throws Exception {
            if (lines < 1 || columns < 1) {
                throw new Exception("Неверный размер.");
            }
            this.lines = lines;
            this.columns = columns;
            this.matrix = new double[lines][columns];
        }

        public Matrix(Matrix init) throws Exception {
            this(init.getLines(), init.getColumns());
            for (int i = 0; i < lines; i++) {
                for (int j = 0; j < columns; j++) {
                    this.matrix[i][j] = init.matrix[i][j];
                }
            }
        }

        public int getLines() {
            return lines;
        }

        public int getColumns() {
            return columns;
        }

        public void setColumn(int index, Vector column) throws Exception {
            if (column.getLength() != lines) {
                throw new Exception("Неверный вектор.");
            }
            for (int i = 0; i < lines; i++) {
                matrix[i][index] = column.vector[i];
            }
        }

        public void print() {
            for (double[] i : matrix) {
                for (double j : i) {
                    System.out.printf("%.5f", j);
                    System.out.print("  ");
                }
                System.out.println();
            }
        }

        public void swap(int fi, int fj, int si, int sj) {
            double tmp = matrix[fi][fj];
            matrix[fi][fj] = matrix[si][sj];
            matrix[si][sj] = tmp;
        }

        public void swapLines(int fline, int sline) {
            for (int i = 0; i < columns; i++) {
                swap(fline, i, sline, i);
            }
        }

        public void fillDefault() {
            double[][] a = {{0.6444, 0.0000, -0.1683, 0.1184, 0.1973},
                    {-0.0395, 0.4208, 0.0000, -0.0802, 0.0263},
                    {0.0132, -0.1184, 0.7627, 0.0145, 0.0460},
                    {0.0395, 0.0000, -0.0960, 0.7627, 0.0000},
                    {0.0263, -0.0395, 0.1907, -0.0158, 0.5523}};
            this.lines = 5;
            this.columns = 5;
            this.matrix = a;
        }

        public Vector mul(Vector vector) throws Exception {
            if (columns != vector.getLength()) {
                throw new Exception("Неверная матрица или вектор.");
            }
            Vector result = new Vector(vector.getLength());
            for (int i = 0; i < lines; i++) {
                result.vector[i] = 0;
                for (int j = 0; j < columns; j++) {
                    result.vector[i] += matrix[i][j] * vector.vector[j];
                }
            }
            return result;
        }

        public Matrix mul(Matrix mtr) throws Exception {
            if (columns != mtr.getLines()) {
                throw new Exception("Неверная матрица.");
            }
            Matrix result = new Matrix(lines, mtr.getColumns());
            for (int i = 0; i < result.getLines(); i++) {
                for (int j = 0; j < result.getColumns(); j++) {
                    result.matrix[i][j] = 0;
                    for (int k = 0; k < columns; k++) {
                        result.matrix[i][j] += this.matrix[i][k] * mtr.matrix[k][j];
                    }
                }
            }
            return result;
        }

        public Matrix transpose() throws Exception {
            if (lines != columns) {
                throw new Exception("Неверная матрица.");
            }
            Matrix result = new Matrix(this);
            for (int i = 0; i < lines; i++) {
                for (int j = i + 1; j < columns; j++) {
                    result.swap(i, j, j, i);
                }
            }
            return result;
        }
    }

    private static class Vector {
        public double[] vector;
        private int length;

        public Vector(int length) throws Exception {
            if (length < 1) {
                throw new Exception("Неверный размер.");
            }
            this.length = length;
            vector = new double[length];
        }

        public int getLength() {
            return length;
        }

        public void print(boolean exponent) {
            for (double item : vector) {
                if (exponent) {
                    System.out.printf("%e\n", item);
                } else {
                    System.out.printf("%.5f\n", item);
                }
            }
        }

        public void swap(int i, int j) {
            double tmp = vector[i];
            vector[i] = vector[j];
            vector[j] = tmp;
        }

        public Vector mul(double num) throws Exception {
            Vector result = new Vector(length);
            for (int i = 0; i < length; i++) {
                result.vector[i] = this.vector[i] * num;
            }
            return result;
        }

        public Vector subtract(Vector sub) throws Exception {
            if (length != sub.getLength()) {
                throw new Exception("Неверный вектор.");
            }
            Vector result = new Vector(length);
            for (int i = 0; i < length; i++) {
                result.vector[i] = this.vector[i] - sub.vector[i];
            }
            return result;
        }

        public double normI() {
            double max = Math.abs(vector[0]);
            for (int i = 1; i < length; i++) {
                if (Math.abs(vector[i]) > max) {
                    max = Math.abs(vector[i]);
                }
            }
            return max;
        }
    }

    private static Matrix A;
    private static Matrix C;
    private static final int n = 5;
    private static double lambda = 0.780861;

    public static void main(String[] args) {
        double[] p;
        Vector eigen;
        Vector r;
        try {
            A = new Matrix(n, n);
            C = new Matrix(n, n);
            A.fillDefault();
            A = A.transpose().mul(A);
            System.out.println("Матрица A:");
            A.print();
            System.out.println();
            p = methodKrylov();
            System.out.println("Коэфициенты: ");
            for (int i = 0; i < n; i++) {
                System.out.print("P" + i + " = ");
                System.out.format("%.5f", p[i]);
                System.out.println();
            }
            System.out.println();
            eigen = findEigenvector(lambda, p);
            System.out.println("Собственный вектор соответствующий максимальному собственному значению " + lambda + " :");
            eigen.print(false);
            System.out.println();
            r = A.mul(eigen).subtract(eigen.mul(lambda));
            System.out.println("Вектор невязки: ");
            r.print(true);
            System.out.println();
            System.out.println("Норма вектора невязки: " + r.normI());
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    private static double[] methodKrylov() throws Exception {
        Vector c = new Vector(n);
        double[] result = new double[n];
        c.vector[0] = 1;
        C.setColumn(n - 1, c);
        for (int i = 1; i < n; i++) {
            c = A.mul(c);
            C.setColumn(n - 1 - i, c);
        }
        return gauss(new Matrix(C), A.mul(c));
    }

    private static double[] gauss(Matrix a, Vector b) throws Exception {
        double[] x = new double[n];
        double max;
        int maxk;
        for (int k = 0; k < n; k++) {
            max = a.matrix[k][k];
            maxk = k;
            for (int i = k + 1; i < n; i++) {
                if (Math.abs(max) < Math.abs(a.matrix[i][k])) {
                    max = a.matrix[i][k];
                    maxk = i;
                }
            }
            if (maxk != k) {
                a.swapLines(k, maxk);
                b.swap(k, maxk);
            }
            for (int j = k; j < n; j++) {
                a.matrix[k][j] /= max;
            }
            b.vector[k] /= max;
            for (int i = k + 1; i < n; i++) {
                for (int j = k + 1; j < n; j++) {
                    a.matrix[i][j] -= a.matrix[i][k] * a.matrix[k][j];
                }
                b.vector[i] -= a.matrix[i][k] * b.vector[k];
                a.matrix[i][k] = 0;
            }
        }
        for (int i = n - 1; i >= 0; i--) {
            x[i] = b.vector[i];
            for (int j = i + 1; j < n; j++) {
                x[i] -= a.matrix[i][j] * x[j];
            }
        }
        return x;
    }

    private static Vector findEigenvector(double lambda, double[] koefs) throws Exception {
        Vector result = new Vector(n);
        Vector beta = new Vector(n);
        for(int i = 0; i < n; i++) {
            beta.vector[i] = Math.pow(lambda, i);
            for(int j = 0; j < i; j++) {
                beta.vector[i] -= koefs[j] * Math.pow(lambda, i - j - 1);
            }
        }
        for(int i = 0; i < n; i++) {
            result.vector[i] = 0;
            for(int j = 0; j < n; j++) {
                result.vector[i] += beta.vector[j] * C.matrix[i][j];
            }
        }
        return result;
    }
}
