package com.unicamp;

import java.math.BigDecimal;
import java.math.RoundingMode;
import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Projeto 2 - MS211 Turma Z
 * Segundo semestre de 2019
 *
 * Integrantes:
 *  - Júlio Moreira Blás de Barros | RA 200491
 *  - Flávio Murilo Reginato       | RA 197088
 */

public class Main {

    public static void main(String[] args) {
        // a)
        var r = 0.5; var K = 10;
        var y_0 = 1; var h = 0.05;

        solveWithEulerAndPrintResults(r, K, y_0, h, 0, 4, true);

        // b) Variando parâmetros:
        r = 0.5; K = 10;
        y_0 = 1; h = 0.01;

        System.out.println("[Decreasing h]");
        solveWithEulerAndPrintResults(r, K, y_0, h, 0, 4, true);
        // Quanto menor o h, maior a precisão

        // Variando parâmetros:
        r = 0.5; K = 100;
        y_0 = 1; h = 0.05;

        System.out.println("[Increasing K]");
        solveWithEulerAndPrintResults(r, K, y_0, h, 0, 4, true);
        // Com a mudança do K, não há mudança muito expressiva no resultado

        // Variando parâmetros:
        r = 0.5; K = 10;
        y_0 = 5; h = 0.05;

        System.out.println("[Increasing y_0]");
        solveWithEulerAndPrintResults(r, K, y_0, h, 0, 4, true);
        // Com o aumento do y_0, houve um aumento de precisão,
        // provavelmente pelo fato dos resultados terem maior valor
        // e, por conta disso, proporcionalmente o erro ser menor.


        // c) Runge-Kutta de quarta ordem
        r = 0.5; K = 10;
        y_0 = 1; h = 0.05;
        solveWithEulerAndPrintResults(r, K, y_0, h, 0, 4, false);
        solveWithKuttaAndPrintResults(r, K, y_0, h, 0, 4, true);
        // Runge kutta de quarta ordem é muito mais preciso, dada a maior complexidade
        // por passo da resolução e a média ponderada


        // Intervalo [0,10]
        r = 0.5; K = 10;
        y_0 = 1; h = 0.05;
        solveWithEulerAndPrintResults(r, K, y_0, h, 0, 10, false);
        solveWithKuttaAndPrintResults(r, K, y_0, h, 0, 10, true);
        // Mesmo em um intervalo maior, o método de Runge Kutta mostrou-se mais preciso,
        // com aproximadamente 6 casas decimais de precisão.

        // d) Definitivamente o método de Runge Kutta é o mais eficiente, com uma precisão de 6 casas decimais
        // Enquanto o método de euler chega a, no máximo, 3.
    }

    private static void solveWithEulerAndPrintResults(double r, int k, int y_0, double h, double start, double end, boolean printAnalytical) {
        var interval = new Interval(start, end);
        var derivative = new VerhulstLogisticDerivative(r, k);
        var solution = new VerhulstLogisticSolution(r, k, y_0);

        System.out.println("Euler method: ");
        System.out.println(EDOSolver.eulerSolve(derivative, interval, h, y_0));
        System.out.println();
        if (printAnalytical) {
            solveAnalyticallyAndPrintResults(h, interval, solution);
        }
    }

    private static void solveWithKuttaAndPrintResults(double r, int k, int y_0, double h, double start, double end, boolean printAnalytical) {
        var interval = new Interval(start, end);
        var derivative = new VerhulstLogisticDerivative(r, k);
        var solution = new VerhulstLogisticSolution(r, k, y_0);

        System.out.println("Runge-Kutta method: ");
        System.out.println(EDOSolver.rungeKuttaSolve(derivative, interval, h, y_0));
        System.out.println();
        if (printAnalytical) {
            solveAnalyticallyAndPrintResults(h, interval, solution);
        }
    }

    private static void solveAnalyticallyAndPrintResults(double h, Interval interval, VerhulstLogisticSolution solution) {
        System.out.println("Analytic Method: ");
        System.out.println(EDOSolver.analyticSolve(solution, interval, h));
        System.out.println();
        System.out.println();
    }
}

class EDOSolver {
    static List<Point> eulerSolve(
        DoubleParamDerivative function, Interval interval,
        /* h: */double step, double y_0
    ) {
        final var range = interval.range(step);
        var lastValue = y_0;
        final var newList = new ArrayList<Point>(range.size());

        for (double x : range) {
            newList.add(new Point(x, MathExtension.round(lastValue)));
                // the first value will be the (x_0, y_0)

            lastValue = lastValue + step * function.evaluate(x, lastValue);
        }

        return newList;
    }

    static List<Point> rungeKuttaSolve(
        DoubleParamDerivative function, Interval interval,
        /* h: */double step, double y_0
    ) {
        final var range = interval.range(step);
        var lastValue = y_0;
        final var newList = new ArrayList<Point>(range.size());

        for (double x : range) {
            newList.add(new Point(x, MathExtension.round(lastValue)));
            // the first value will be the (x_0, y_0)

            final var k1 = function.evaluate(x, lastValue);
            final var k2 = function.evaluate(x + step/2, lastValue + k1 * step/2);
            final var k3 = function.evaluate(x + step/2, lastValue + k2 * step/2);
            final var k4 = function.evaluate(x + step, lastValue + k3 * step);
            final var k = (1.0/6.0) * (k1 + 2*k2 + 2*k3 + k4);

            lastValue = lastValue + step * k;
        }

        return newList;
    }

    static List<Point> analyticSolve(
        SingleParamFunction function, Interval interval,/* h: */double step
    ) {
        final var range = interval.range(step);

        return range.stream()
            .map(x -> new Point(x, MathExtension.round(function.evaluate(x))))
            .collect(Collectors.toList());
    }
}


interface SingleParamFunction {
    double evaluate(double x);
}

interface DoubleParamDerivative {
    double evaluate(double x, double y);
}

class VerhulstLogisticDerivative implements DoubleParamDerivative {

    private final double r;
    private final double K;

    VerhulstLogisticDerivative(double r, double k) {
        this.r = r;
        this.K = k;
    }

    @Override
    public double evaluate(double x, double y) {
        return r * y * (1 - (y/K));
    }
}

class VerhulstLogisticSolution implements SingleParamFunction {

    private final double r;
    private final double K;
    private final double y_0;

    VerhulstLogisticSolution(double r, double k, double y_0) {
        this.r = r;
        this.K = k;
        this.y_0 = y_0;
    }

    @Override
    public double evaluate(double x) {
        var exp = Math.exp(r * x);
        return (K * y_0 * exp) / (K + y_0*(exp - 1));
    }
}

// Data classes:

class Interval {
    private static final double EPSILON = 0.000001;
    private final double start;
    private final double end;

    Interval(double start, double end) {
        this.start = start;
        this.end = end;
    }

    public double getStart() {
        return start;
    }

    public double getEnd() {
        return end;
    }

    List<Double> range(double step) {
        final int maxSize = (int) Math.round(Math.ceil( (end-start)/step ));
        final var range = new ArrayList<Double>(maxSize);

        for (double i = start; (i - end) <= EPSILON; i += step) {
            range.add(MathExtension.round(i));
        }

        return range;
    }
}

class Point {
    private final double x;
    private final double y;

    Point(double x, double y) {
        this.x = x;
        this.y = y;
    }

    public double getX() {
        return x;
    }

    public double getY() {
        return y;
    }

    public String toString() {
        return "("+x+", "+y+")";
    }
}

class MathExtension {
    public static final int GENERAL_PRECISION = 10;
    public static double round(double value) {
        BigDecimal number = BigDecimal.valueOf(value).setScale(GENERAL_PRECISION, RoundingMode.HALF_UP);
        return number.doubleValue();
    }
}