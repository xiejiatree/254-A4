/*
    MST.java

    Euclidean minimum spanning tree implementation.

    Uses Dwyer's algorithm for Delaunay triangularization, then
    Kruskal's MST algorithm on the resulting mesh.  You need to
    parallelize both stages.  For the triangulation stage, I recommend
    creating a new thread class similar to worker, to use in the
    divide-and-conquer step of triangulate().  For the tree stage, I
    recommend letting worker threads find subtrees to merge in parallel,
    but forcing them to finalize the merges in order (after
    double-checking to make sure their subtrees haven't been changed by
    earlier finalizations).

    There are better ways to parallelize each of these steps, but
    they're quite a bit harder.

     Accepts the following command-line arguments:

    −a  [0123]
        Animation mode. 
        0   (default) =>
            print run time to standard output, but nothing else 
        1 =>
            print list of created, destroyed, and selected (tree) edges,
            plus run time 
        2 =>
            create a GUI that shows the triangulation and MST, and allow
            the user to re-run with additional sets of points 
        3 =>
            animate the algorithm on the screen as it runs.  
    −n  num
        Number of points.  Default = 50.  More than a couple hundred becomes
        too dense to look good when animated.  You’ll need to run big numbers
        (more than 10,000) to get multi-second execution times.  
    −s  num
        Seed for pseudorandom number generator.  Every value of the seed
        produces a different set of points.  
    −t  num
        Number of threads (max) that should be running at any given time.
        This argument is currently unused; it’s it’s here to support your
        parallelization efforts.  

    Code currently assumes (without checking) that no three points are
    co-linear and no four points are co-circular.

    (c) Michael L. Scott, 2025.  Based heavily on earlier incarnations
    of several programming projects, and on Delaunay mesh code written
    in 2007.  For use by students in CSC 2/454 at the University of
    Rochester, during the spring 2025 term.  All other use requires
    written permission of the author.
 */

import java.awt.*;
import java.awt.event.*;
import java.io.*;
import javax.swing.*;
import java.util.*;
import java.lang.*;
import java.util.concurrent.ConcurrentSkipListSet;

public class MST {
    private static int n = 50;              // default number of points
    private static long sd = 0;             // default random number seed
    private static int numThreads = 1;      // default
        // This variable is not yet used; you'll want to use it for your
        // parallel version.

    private static final int TIMING_ONLY    = 0;
    private static final int PRINT_EVENTS   = 1;
    private static final int SHOW_RESULT    = 2;
    private static final int FULL_ANIMATION = 3;
    private static int animationMode = TIMING_ONLY;       // default

    // Process command-line arguments.
    //
    private static void parseArgs(String[] args) {
        for (int i = 0; i < args.length; i++) {
            if (args[i].equals("-a")) {
                if (++i >= args.length) {
                    System.err.print("Missing animation level\n");
                } else {
                    int an = -1;
                    try {
                        an = Integer.parseInt(args[i]);
                    } catch (NumberFormatException e) { }
                    if (an >= TIMING_ONLY && an <= FULL_ANIMATION) {
                        animationMode = an;
                    } else {
                        System.err.printf("Invalid animation level: %s\n", args[i]);
                    }
                }
            } else if (args[i].equals("-n")) {
                if (++i >= args.length) {
                    System.err.print("Missing number of points\n");
                } else {
                    int np = -1;
                    try {
                        np = Integer.parseInt(args[i]);
                    } catch (NumberFormatException e) { }
                    if (np > 0) {
                        n = np;
                    } else {
                        System.err.printf("Invalid number of points: %s\n", args[i]);
                    }
                }
            } else if (args[i].equals("-s")) {
                if (++i >= args.length) {
                    System.err.print("Missing seed\n");
                } else {
                    try {
                        sd = Integer.parseInt(args[i]);
                    } catch (NumberFormatException e) {
                        System.err.printf("Invalid seed: %s\n", args[i]);
                    }
                }
            } else if (args[i].equals("-t")) {
                if (++i >= args.length) {
                    System.err.print("Missing number of threads\n");
                } else {
                    int nt = -1;
                    try {
                        nt = Integer.parseInt(args[i]);
                    } catch (NumberFormatException e) { }
                    if (nt > 0) {
                        numThreads = nt;
                    } else {
                        System.err.printf("Invalid number of threads: %s\n", args[i]);
                    }
                }
            } else {
                System.err.printf("Unexpected argument: %s\n", args[i]);
            }
        }
    }

    // Initialize program components for specified animation mode.
    //
    private MSTworld build(RootPaneContainer pane, int an) {
        final Coordinator c = new Coordinator();
        MSTworld w = new MSTworld(n, sd, c);
        Animation t = null;
        if (an == SHOW_RESULT || an == FULL_ANIMATION) {
            t = new Animation(w);
            new UI(c, w, t, sd, pane);
        }
        final Animation a = t;
        if (an == PRINT_EVENTS) {
            w.setHooks(
                new MSTworld.EdgeRoutine() {
                    public void run(int x1, int y1, int x2, int y2, boolean dum) {
                        System.out.printf("created   %12d %12d %12d %12d\n",
                                          x1, y1, x2, y2);
                    }},
                new MSTworld.EdgeRoutine() {
                    public void run(int x1, int y1, int x2, int y2, boolean dum) {
                        System.out.printf("destroyed %12d %12d %12d %12d\n",
                                          x1, y1, x2, y2);
                    }},
                new MSTworld.EdgeRoutine() {
                    public void run(int x1, int y1, int x2, int y2, boolean dum) {
                        System.out.printf("selected  %12d %12d %12d %12d\n",
                                          x1, y1, x2, y2);
                    }});
        } else if (an == FULL_ANIMATION) {
            MSTworld.EdgeRoutine er = new MSTworld.EdgeRoutine() {
                public void run(int x1, int y1, int x2, int y2, boolean dum)
                        throws Coordinator.KilledException {
                    c.hesitate();
                    a.repaint();        // graphics need to be re-rendered
                }};
            w.setHooks(er, er, er);
        }
        return w;
    }

    public static void main(String[] args) {
        parseArgs(args);
        MST me = new MST();
        JFrame f = null;
        if (animationMode == SHOW_RESULT || animationMode == FULL_ANIMATION) {
            f = new JFrame("MST");
            f.addWindowListener(new WindowAdapter() {
                public void windowClosing(WindowEvent e) {
                    System.exit(0);
                }
            });
        } else {
            System.out.printf("%d points, seed %d\n", n, sd);
        }
        MSTworld w = me.build(f, animationMode);
        if (f != null) {
            f.pack();
            f.setVisible(true);
        } else {
            // Using terminal I/O rather than graphics.
            // Execute the guts of the run button handler method here.
            long startTime = new Date().getTime();
            long midTime = 0;
            try {
                w.DwyerSolve();
                midTime = new Date().getTime();
                w.KruskalSolve();
            } catch(Coordinator.KilledException e) { }
            long endTime = new Date().getTime();
            System.out.printf("elapsed time: %.3f + %.3f = %.3f seconds\n",
                              (double) (midTime-startTime)/1000,
                              (double) (endTime-midTime)/1000,
                              (double) (endTime-startTime)/1000);
        }
    }
}

// The Worker is the thread that does the actual work of finding a
// triangulation and MST (in the animated version -- main thread does it
// in the terminal I/O version).
//
class Worker extends Thread {
    private final MSTworld w;
    private final Coordinator c;
    private final UI u;
    private final Animation a;

    // Constructor
    //
    public Worker(MSTworld w, Coordinator c, UI u, Animation a) {
        this.w = w;
        this.c = c;
        this.u = u;
        this.a = a;
    }

    // The run() method of a Java Thread is never invoked directly by
    // user code. Rather, it is called by the Java runtime when user
    // code calls start().
    //
    // The run() method of a worker thread *must* begin by calling
    // c.register() and end by calling c.unregister().  These allow the
    // user interface (via the Coordinator) to pause and terminate
    // workers.  Note how the worker is set up to catch KilledException.
    // In the process of unwinding back to here we'll cleanly and
    // automatically release any monitor locks.  If you create new kinds
    // of workers (as part of a parallel solver), make sure they call
    // c.register() and c.unregister() properly.
    //
    public void run() {
        try {
            c.register();
            w.DwyerSolve();
            w.KruskalSolve();
            c.unregister();
        } catch(Coordinator.KilledException e) { }
        if (a != null) {
            // Tell the graphics event thread to unset the default
            // button when it gets a chance.  (Threads other than the
            // event thread cannot safely modify the GUI directly.)
            a.repaint();
            SwingUtilities.invokeLater(new Runnable() {
                public void run() {
                    u.getRootPane().setDefaultButton(null);
                    u.updateTime();
                }
            });
        }
    }
}

// The MSTworld contains all points and edges.
//
class MSTworld {
    // Much of the logic in Dwyer's algorithm is parameterized by
    // directional (X or Y) and rotational (clockwise, counterclockwise)
    // orientation.  The following constants get plugged into the
    // parameter slots.
    private static final int xdim = 0;
    private static final int ydim = 1;
    private static final int ccw = 0;
    private static final int cw = 1;

    private int minx;   // smallest x value among all points
    private int maxx;   // largest x value among all points
    private int miny;   // smallest y value among all points
    private int maxy;   // largest y value among all points
    public int getMinx() {return minx;}
    public int getMaxx() {return maxx;}
    public int getMiny() {return miny;}
    public int getMaxy() {return maxy;}

    // The following 7 fields are set by the MSTworld constructor.
    private final Coordinator coord;
        // Not needed at present, but will need to be passed to any
        // newly created workers.
    private final int n;  // number of points
    private final point points[];
        // main array of points, used for partitioning and rendering
    private final HashSet<point> pointHash;
        // Used to ensure that we never have two points directly on top of
        // each other.  See point.hashCode and point.equals below.
    private final SortedSet<edge> edges;
        // Used for rendering.  Ordering supports the KruskalSolve stage.
    private long sd = 0;
    private final Random prn;     // pseudo-random number generator

    // Constructor
    //
    public MSTworld(int n, long sd, Coordinator c) {
        this.n = n;
        this.sd = sd;
        coord = c;

        points = new point[n];
        edges = new ConcurrentSkipListSet<edge>(new edgeComp());
            // Supports safe concurrent access by worker and graphics threads,
            // and as a SortedSet it keeps the edges in order by length.
        pointHash = new HashSet<point>(n);

        prn = new Random();
        reset();
    }

    // 3x3 determinant.  Called by 4x4 determinant.
    //
    private double det3(double a, double b, double c,
                        double d, double e, double f,
                        double g, double h, double i) {
        return a * (e*i - f*h)
             - b * (d*i - f*g)
             + c * (d*h - e*g);
    }

    // 4x4 determinant.  Called by encircled (below).
    //
    private double det4(double a, double b, double c, double d,
                        double e, double f, double g, double h,
                        double i, double j, double k, double l,
                        double m, double n, double o, double p) {
        return a * det3(f, g, h, j, k, l, n, o, p)
             - b * det3(e, g, h, i, k, l, m, o, p)
             + c * det3(e, f, h, i, j, l, m, n, p)
             - d * det3(e, f, g, i, j, k, m, n, o);
    }

    // If A, B, and C are on a circle, in counter-clockwise order, then
    // D lies within that circle iff the following determinant is positive:
    //
    // | Ax  Ay  Ax^2+Ay^2  1 |
    // | Bx  By  Bx^2+By^2  1 |
    // | Cx  Cy  Cx^2+Cy^2  1 |
    // | Dx  Dy  Dx^2+Dy^2  1 |
    //
    private boolean encircled(point A, point B, point C, point D, int dir) {
        if (dir == cw) {
            point t = A;  A = C;  C = t;
        }
        double Ax = A.getCoord(xdim);   double Ay = A.getCoord(ydim);
        double Bx = B.getCoord(xdim);   double By = B.getCoord(ydim);
        double Cx = C.getCoord(xdim);   double Cy = C.getCoord(ydim);
        double Dx = D.getCoord(xdim);   double Dy = D.getCoord(ydim);

        return det4(Ax, Ay, (Ax*Ax + Ay*Ay), 1,
                    Bx, By, (Bx*Bx + By*By), 1,
                    Cx, Cy, (Cx*Cx + Cy*Cy), 1,
                    Dx, Dy, (Dx*Dx + Dy*Dy), 1) > 0;
    }

    // swap points[i] and points[j]
    //
    private void swap(int i, int j) {
        point t = points[i];
        points[i] = points[j];
        points[j] = t;
    }

    // A point is a mesh/tree vertex.  It also serves, in the Kruskal
    // stage, as a union-find set.  Its x and y coordinates are private
    // and final.  Use the getCoord method to read their values.
    //
    private class point {
        private final int coordinates[] = new int[2];
        public edge firstEdge;

        // The following two fields are needed in the Kruskal stage (only).
        private point representative = null;    // equivalence set (subtree)
        private int subtreeSize = 1;

        // Constructor
        //
        public point(int x, int y) {
            coordinates[xdim] = x;  coordinates[ydim] = y;
            // firstEdge == null
        }

        // This is essentially the "find" of a union-find implementation.
        public point subtree() {
            point p = this;
            while (p.representative != null) p = p.representative;
            return p;
        }

        // And this is the "union".
        public void merge(point p) {
            // Make larger set the representative of the smaller set.
            assert representative == null;
            if (subtreeSize > p.subtreeSize) {
                subtreeSize += p.subtreeSize;
                p.representative = this;
            } else {
                p.subtreeSize += subtreeSize;
                representative = p;
            }
        }

        private int getCoord(int dim) {
            return coordinates[dim];
        }

        // Override Object.hashCode and Object.equals.  This way two
        // points are equal (and hash to the same slot in HashSet
        // pointHash) if they have the same coordinates, even if they
        // are different objects.
        //
        public int hashCode() {
            return coordinates[xdim] ^ coordinates[ydim];
        }
        public boolean equals(Object o) {
            point p = (point) o;            // run-time type check
            return p.coordinates[xdim] == coordinates[xdim]
                && p.coordinates[ydim] == coordinates[ydim];
        }
    }

    // Signatures for things someone might want us to do with a point or
    // an edge (e.g., display it).
    // 
    public interface PointRoutine{
        public void run(int x, int y);
    }
    public interface EdgeRoutine {
        public void run(int x1, int y1, int x2, int y2, boolean treeEdge)
            throws Coordinator.KilledException;
    }

    public void forAllPoints(PointRoutine pr) {
        for (point p : points) {
            pr.run(p.getCoord(xdim), p.getCoord(ydim));
        }
    }
    public void forAllEdges(EdgeRoutine er) {
        for (edge e : edges) {
            try {
                er.run(e.points[0].getCoord(xdim),
                       e.points[0].getCoord(ydim),
                       e.points[1].getCoord(xdim),
                       e.points[1].getCoord(ydim), e.isMSTedge);
            } catch (Coordinator.KilledException f) { }
        }
    }

    // Routines to call when performing the specified operations:
    private static EdgeRoutine edgeCreateHook = null;
    private static EdgeRoutine edgeSelectHook = null;
    private static EdgeRoutine edgeDestroyHook = null;

    // The following is separate from the constructor to avoid a
    // circularity problem: when working in FULL_ANIMATION mode, the
    // Animation object needs a reference to the MSTworld object, and the
    // MSTworld object needs references to the hooks of the Animation object.
    //
    public void setHooks(EdgeRoutine ech, EdgeRoutine edh, EdgeRoutine esh) {
        edgeCreateHook = ech;
        edgeDestroyHook = edh;
        edgeSelectHook = esh;
    }

    // Edges encapsulate the bulk of the information about the triangulation.
    // Each edge contains references to its endpoints and to the next
    // edges clockwise and counterclockwise about those endpoints.
    //
    private class edge {
        public final point[] points = new point[2];
        public final edge[][] neighbors = new edge[2][2];
            // indexed first by edge end and then by rotational direction
        private boolean isMSTedge = false;
        public final double length;

        // Return index of point p within edge
        //
        public int indexOf(point p) {
            if (points[0] == p) return 0;
            if (points[1] == p) return 1;
            return -1;      // so I get an error if I use it
        }

        // utility routine for constructor
        //
        private void initializeEnd(point p, edge e, int end, int dir) {
            if (e == null) {
                neighbors[end][dir] = neighbors[end][1-dir] = this;
                p.firstEdge = this;
            } else {
                int i = e.indexOf(p);
                neighbors[end][1-dir] = e;
                neighbors[end][dir] = e.neighbors[i][dir];
                e.neighbors[i][dir] = this;
                i = neighbors[end][dir].indexOf(p);
                neighbors[end][dir].neighbors[i][1-dir] = this;
            }
        }

        // Constructor: connect points a and b, inserting dir (CW or CCW)
        // of edge ea at the a end and 1-dir of edge eb at the b end.
        // Either or both of ea and eb may be null.
        //
        public edge(point a, point b, edge ea, edge eb, int dir)
                throws Coordinator.KilledException {
            points[0] = a;  points[1] = b;
            double dx = (double) a.getCoord(xdim) - (double) b.getCoord(xdim);
            double dy = (double) a.getCoord(ydim) - (double) b.getCoord(ydim);
            length = Math.sqrt(dx * dx + dy * dy);

            initializeEnd(a, ea, 0, dir);
            initializeEnd(b, eb, 1, 1-dir);

            edges.add(this);
            if (edgeCreateHook != null)
                edgeCreateHook.run(points[0].getCoord(xdim),
                                   points[0].getCoord(ydim),
                                   points[1].getCoord(xdim),
                                   points[1].getCoord(ydim), false);
        }

        // Destructor: take self out of edges, point edge lists.
        // Should only be called when flipping an edge, so destroyed
        // edge should have neighbors at both ends.
        //
        public void destroy() throws Coordinator.KilledException {
            edges.remove(this);
            for (int i = 0; i < 2; i++) {
                int cw_index = neighbors[i][cw].indexOf(points[i]);
                int ccw_index = neighbors[i][ccw].indexOf(points[i]);
                neighbors[i][cw].neighbors[cw_index][ccw] = neighbors[i][ccw];
                neighbors[i][ccw].neighbors[ccw_index][cw] = neighbors[i][cw];
                if (points[i].firstEdge == this)
                    points[i].firstEdge = neighbors[i][ccw];
            }
            if (edgeDestroyHook != null)
                edgeDestroyHook.run(points[0].getCoord(xdim),
                                    points[0].getCoord(ydim),
                                    points[1].getCoord(xdim),
                                    points[1].getCoord(ydim), false);
        }

        // Assume edges are unique.  Override Object.equals to make it
        // consistent with edgeComp.compare below.
        //
        public boolean equals(Object o) {
            return this == o;
        }

        // Label this edge as an MST edge.
        //
        public void addToMST() throws Coordinator.KilledException {
            isMSTedge = true;
            if (edgeSelectHook != null)
                edgeSelectHook.run(points[0].getCoord(xdim),
                                   points[0].getCoord(ydim),
                                   points[1].getCoord(xdim),
                                   points[1].getCoord(ydim), false);
        }
    }

    // To support ordered set of edges.  Return 0 _only_ if two
    // arguments are the _same_ edge (this is necessary for unique
    // membership in set).  Otherwise order based on length.
    // If lengths are the same, order by coordinates.
    //
    public static class edgeComp implements Comparator<edge> {
        public int compare(edge e1, edge e2) {
            if (e1.equals(e2)) return 0;
            if (e1.length < e2.length) return -1;
            if (e1.length > e2.length) return 1;
            int e1xmin = e1.points[0].getCoord(xdim)
                            < e1.points[1].getCoord(xdim) ?
                                e1.points[0].getCoord(xdim) :
                                e1.points[1].getCoord(xdim);
            int e2xmin = e2.points[0].getCoord(xdim)
                            < e2.points[1].getCoord(xdim) ?
                                e2.points[0].getCoord(xdim) :
                                e2.points[1].getCoord(xdim);
            if (e1xmin < e2xmin) return -1;
            if (e1xmin > e2xmin) return 1;
            int e1ymin = e1.points[0].getCoord(ydim)
                            < e1.points[1].getCoord(ydim) ?
                                e1.points[0].getCoord(ydim) :
                                e1.points[1].getCoord(ydim);
            int e2ymin = e2.points[0].getCoord(ydim)
                            < e2.points[1].getCoord(ydim) ?
                                e2.points[0].getCoord(ydim) :
                                e2.points[1].getCoord(ydim);
            if (e1ymin < e2ymin) return -1;
            // if (e1ymin > e2ymin)
            return 1;
            // no other options; endpoints have to be distinct
        }
    }

    // Called by the UI when it wants to reset with a new seed.
    //
    public long randomize() {
        sd++;
        reset();
        return sd;
    }

    // Called by the UI when it wants to start over.
    //
    public void reset() {
        prn.setSeed(sd);
        minx = Integer.MAX_VALUE;   miny = Integer.MAX_VALUE;
        maxx = Integer.MIN_VALUE;   maxy = Integer.MIN_VALUE;
        pointHash.clear();      // empty out the set of points
        for (int i = 0; i < n; i++) {
            point p;
            int x;
            int y;
            do {
                x = prn.nextInt();
                y = prn.nextInt();
                p = new point(x, y);
            } while (pointHash.contains(p));
            pointHash.add(p);
            if (x < minx) minx = x;
            if (y < miny) miny = y;
            if (x > maxx) maxx = x;
            if (y > maxy) maxy = y;
            points[i] = p;
        }
        edges.clear();      // empty out the set of edges
    }

    // Is angle from p1 to p2 to p3, in direction dir
    // around p2, greater than or equal to 180 degrees?
    //
    private boolean externAngle(point p1, point p2, point p3, int dir) {
        if (dir == cw) {
            point t = p1;  p1 = p3;  p3 = t;
        }
        int x1 = p1.getCoord(xdim);     int y1 = p1.getCoord(ydim);
        int x2 = p2.getCoord(xdim);     int y2 = p2.getCoord(ydim);
        int x3 = p3.getCoord(xdim);     int y3 = p3.getCoord(ydim);

        if (x1 == x2) {                         // first segment vertical
            if (y1 > y2) {                      // points down
                return (x3 >= x2);
            } else {                            // points up
                return (x3 <= x2);
            }
        } else {
            double m = (((double) y2) - y1) / (((double) x2) - x1);
                // slope of first segment
            if (x1 > x2) {                      // points left
                return (y3 <= m * (((double) x3) - x1) + y1);
                // p3 below line
            } else {                            // points right
                return (y3 >= m * (((double) x3) - x1) + y1);
                // p3 above line
            }
        }
    }

    // Divide points[l..r] into two partitions.  Solve recursively, then
    // stitch back together.  Dim0 values range from [low0..high0].
    // Dim1 values range from [low1..high1].  We partition based on dim0.
    // Base case when 1, 2, or 3 points.
    //
    // As suggested by Dwyer, we swap axes and rotational directions
    // at successive levels of recursion, to minimize the number of long
    // edges that are likely to be broken when stitching.
    //
    private void triangulate(int l, int r, int low0, int high0,
                             int low1, int high1, int parity)
        throws Coordinator.KilledException {

        final int dim0;  final int dim1;
        final int dir0;  final int dir1;

        if (parity == 0) {
            dim0 = xdim;  dim1 = ydim;
            dir0 = ccw;   dir1 = cw;
        } else {
            dim0 = ydim;  dim1 = xdim;
            dir0 = cw;    dir1 = ccw;
        }

        if (l == r) {
            return;
        }
        if (l == r-1) {
            new edge(points[l], points[r], null, null, dir1);
                // direction doesn't matter in this case
            return;
        }
        if (l == r-2) {     // make single triangle
            edge e2 = new edge(points[l+1], points[r], null, null, dir1);
            edge e1 = new edge(points[l], points[l+1], null, e2, dir1);
            if (externAngle(points[l], points[l+1], points[r], dir0)) {
                // new edge is dir0 of edge 1, dir1 of edge 2
                new edge(points[l], points[r], e1, e2, dir0);
            } else {
                // new edge is dir1 of edge 1, dir0 of edge 2
                new edge(points[l], points[r], e1, e2, dir1);
            }
            return;
        }

        // At this point we know we're not a base case; have to subdivide.

        int mid = low0/2 + high0/2;
        int i = l;  int j = r;

        point lp = points[l];          // rightmost point in left half;
        int lp0 = Integer.MIN_VALUE;   // X coord of lp
        point rp = points[r];          // leftmost point in right half;
        int rp0 = Integer.MAX_VALUE;   // X coord of rp

        while (true) {
            // invariants: [i..j] are unexamined;
            // [l..i) are all <= mid; (j..r] are all > mid.

            int i0 = 0;  int j0 = 0;

            while (i < j) {
                i0 = points[i].getCoord(dim0);
                if (i0 > mid) {     // belongs in right half
                    if (i0 < rp0) {
                        rp0 = i0;  rp = points[i];
                    }
                    break;
                } else {
                    if (i0 > lp0) {
                        lp0 = i0;  lp = points[i];
                    }
                }
                i++;
            }

            while (i < j) {
                j0 = points[j].getCoord(dim0);
                if (j0 <= mid) {    // belongs in left half
                    if (j0 > lp0) {
                        lp0 = j0;  lp = points[j];
                    }
                    break;
                } else {
                    if (j0 < rp0) {
                        rp0 = j0;  rp = points[j];
                    }
                }
                j--;
            }

            // at this point either i == j == only unexamined element
            // or i < j (found elements that need to be swapped)
            // or i = j+1 (and all elements are in order)
            if (i == j) {
                i0 = points[i].getCoord(dim0);
                if (i0 > mid) {
                    // give border element to right half
                    if (i0 < rp0) {
                        rp0 = i0;  rp = points[i];
                    }
                    i--;
                } else {
                    // give border element to left half
                    if (i0 > lp0) {
                        lp0 = i0;  lp = points[i];
                    }
                    j++;
                }
                break;
            }
            if (i > j) {
                i--;  j++;  break;
            }
            swap(i, j);
            i++;  j--;
        }
        // Now [l..i] is the left partition and [j..r] is the right.
        // Either may be empty.

        if (i < l) {
            // empty left half
            triangulate(j, r, low1, high1, mid, high0, 1-parity);
        } else if (j > r) {
            // empty right half
            triangulate(l, i, low1, high1, low0, mid, 1-parity);
        } else {
            // divide and conquer
            triangulate(j, r, low1, high1, mid, high0, 1-parity);
            triangulate(l, i, low1, high1, low0, mid, 1-parity);

            // prepare to stitch meshes together up the middle:
            class side {
                public point p;     // working point
                public edge a;      // above p
                public edge b;      // below p
                public point ap;    // at far end of a
                public point bp;    // at far end of b
                public int ai;      // index of p within a
                public int bi;      // index of p within b
            };
            side right = new side();
            side left = new side();
            right.p = rp;
            left.p = lp;

            // Rotate around extreme point to find edges adjacent to Y
            // axis.  This class is basically a hack to get around the
            // lack of nested subroutines in Java.  We invoke its run
            // method twice below.
            //
            class rotateClass {
                void run(side s, int dir) {
                    // rotate around s.p to find edges adjacent to Y axis
                    if (s.p.firstEdge != null) {
                        s.a = s.p.firstEdge;
                        s.ai = s.a.indexOf(s.p);
                        s.ap = s.a.points[1-s.ai];
                        if (s.a.neighbors[s.ai][dir] == s.a) {
                            // only one incident edge on the right
                            s.b = s.a;
                            s.bi = s.ai;
                            s.bp = s.ap;
                        } else {
                            // >= 2 incident edges on the right;
                            // need to find correct ones
                            while (true) {
                                s.b = s.a.neighbors[s.ai][dir];
                                s.bi = s.b.indexOf(s.p);
                                s.bp = s.b.points[1-s.bi];
                                if (externAngle(s.ap, s.p, s.bp, dir)) break;
                                s.a = s.b;
                                s.ai = s.bi;
                                s.ap = s.bp;
                            }
                        }
                    }
                }
            }
            rotateClass rotate = new rotateClass();
            rotate.run(right, dir0);
            rotate.run(left, dir1);

            // Find endpoint of bottom edge of seam, by moving around border
            // as far as possible without going around a corner.  This, too,
            // is basically a nested subroutine.
            //
            class findBottomClass {
                boolean move(side s, int dir, point o) {
                    boolean progress = false;
                    if (s.b != null) {
                        while (!externAngle(s.bp, s.p, o, 1-dir)) {
                            // move s.p in direction dir
                            progress = true;
                            s.a = s.b;
                            s.ai = 1-s.bi;
                            s.ap = s.p;
                            s.p = s.b.points[1-s.bi];
                            s.b = s.b.neighbors[1-s.bi][dir];
                            s.bi = s.b.indexOf(s.p);
                            s.bp = s.b.points[1-s.bi];
                        }
                    }
                    return progress;
                }
            }
            findBottomClass findBottom = new findBottomClass();
            do {} while (findBottom.move(left, dir1, right.p)
                      || findBottom.move(right, dir0, left.p));

            // create bottom edge:
            edge base = new edge(left.p, right.p,
                                 left.a == null ? left.b : left.a,
                                 right.a == null ? right.b : right.a,
                                 dir1);
            final edge bottom = base;
            if (left.a == null) left.a = bottom;
                // left region is a singleton
            if (right.a == null) right.a = bottom;
                // right region is a singleton

            // Work up the seam creating new edges and deleting old
            // edges where necessary.  Note that {left,right}.{b,bi,bp}
            // are no longer needed.

            while (true) {

                // Find candidate endpoint.  Yet another nested subroutine.
                //
                class findCandidateClass {
                    point call(side s, int dir, edge base, point o)
                            throws Coordinator.KilledException {
                            // o is at far end of base
                        if (s.a == bottom) {
                            // region is a singleton
                            return null;
                        }
                        point c = s.a.points[1-s.ai];
                        if (externAngle(o, s.p, c, dir)) {
                            // no more candidates
                            return null;
                        }
                        while (true) {
                            edge na = s.a.neighbors[s.ai][dir];
                                // next edge into region
                            if (na == base) {
                                // wrapped all the way around
                                return c;
                            }
                            int nai = na.indexOf(s.p);
                            point nc = na.points[1-nai];
                                // next potential candidate
                            if (encircled(o, c, s.p, nc, dir)) {
                                // have to break an edge
                                s.a.destroy();
                                s.a = na;
                                s.ai = nai;
                                c = nc;
                            } else return c;
                        }
                    }
                }
                findCandidateClass findCandidate = new findCandidateClass();
                point lc = findCandidate.call(left, dir0, bottom, right.p);
                point rc = findCandidate.call(right, dir1, bottom, left.p);

                if (lc == null && rc == null) {
                    // no more candidates
                    break;
                }
                // Choose between candidates:
                if (lc != null && rc != null &&
                        encircled (right.p, lc, left.p, rc, dir0)) {
                    // Left candidate won't work; circumcircle contains
                    // right candidate.
                    lc = null;
                }
                // Now we know one candidate is null and the other is not.
                if (lc == null) {
                    // use right candidate
                    right.a = right.a.neighbors[1-right.ai][dir1];
                    right.ai = right.a.indexOf(rc);
                    right.ap = right.a.points[1-right.ai];
                    right.p = rc;
                    base = new edge(left.p, rc, left.a, right.a, dir1);
                } else {
                    // use left candidate
                    left.a = left.a.neighbors[1-left.ai][dir0];
                    left.ai = left.a.indexOf(lc);
                    left.ap = left.a.points[1-left.ai];
                    left.p = lc;
                    base = new edge(lc, right.p, left.a, right.a, dir1);
                }
            }
        }
    }

    // This is the actual MST calculation.
    //
    public void KruskalSolve()
        throws Coordinator.KilledException {
        int numTrees = n;
        for (edge e : edges) {
            point st1 = e.points[0].subtree();
            point st2 = e.points[1].subtree();
            if (st1 != st2) {
                // This edge joins two previously separate subtrees.
                st1.merge(st2);
                e.addToMST();
                if (--numTrees == 1) break;
            }
        }
    }

    // This is a wrapper for the root call to triangulate().
    //
    public void DwyerSolve() throws Coordinator.KilledException {
        triangulate(0, n-1, minx, maxx, miny, maxy, 0);
    }
}

// Class Animation is the one really complicated sub-pane of the user interface.
// 
class Animation extends JPanel {
    private static final int width = 512;      // canvas dimensions
    private static final int height = 512;
    private static final int dotsize = 6;
    private static final int border = dotsize;
    private final MSTworld w;

    // Constructor
    //
    public Animation(MSTworld w) {
        setPreferredSize(new Dimension(width+border*2, height+border*2));
        setBackground(Color.white);
        setForeground(Color.black);
        this.w = w;
        reset();
    }

    // The next two routines figure out where to render the dot
    // for a point, given the size of the animation panel and the spread
    // of x and y values among all points.
    //
    private int xPosition(int x) {
        return (int)
            (((double)x - (double)w.getMinx()) * (double)width
                / ((double)w.getMaxx() - (double)w.getMinx())) + border;
    }
    private int yPosition(int y) {
        return (int)
            (((double)w.getMaxy() - (double)y) * (double)height
                / ((double)w.getMaxy() - (double)w.getMiny())) + border;
    }

    // The following method is called automatically by the graphics
    // system when it thinks the Animation canvas needs to be
    // re-displayed.  This can happen because code elsewhere in this
    // program called repaint(), or because of hiding/revealing or
    // open/close operations in the surrounding window system.
    //
    public void paintComponent(final Graphics g) {
        final Graphics2D g2 = (Graphics2D) g;

        super.paintComponent(g);    // clears panel
        w.forAllEdges(new MSTworld.EdgeRoutine() {
            public void run(int x1, int y1, int x2, int y2, boolean bold) {
                if (bold) {
                    g2.setPaint(Color.red);
                    g2.setStroke(new BasicStroke(4));
                } else {
                    g2.setPaint(Color.gray);
                    g2.setStroke(new BasicStroke(1));
                }
                g.drawLine(xPosition(x1), yPosition(y1),
                           xPosition(x2), yPosition(y2));
            }
        });
        w.forAllPoints(new MSTworld.PointRoutine() {
            public void run(int x, int y) {
                g2.setPaint(Color.blue);
                g.fillOval(xPosition(x)-dotsize/2, yPosition(y)-dotsize/2,
                           dotsize, dotsize);
            }
        });
    }

    // UI needs to call this routine when point locations have changed.
    //
    public void reset() {
        repaint();      // Tell graphics system to re-render.
    }
}

// Class UI is the user interface.  It displays an MSTworld canvas above
// a row of buttons and a row of statistics.  Actions (event handlers)
// are defined for each of the buttons.  Depending on the state of the
// UI, either the "run" or the "pause" button is the default (highlighted in
// most window systems); it will often self-push if you hit carriage return.
//
class UI extends JPanel {
    private final Coordinator coordinator;
    private final MSTworld world;
    private final Animation animation;

    private final JRootPane root;
    private static final int externalBorder = 6;

    private static final int stopped = 0;
    private static final int running = 1;
    private static final int paused = 2;

    private int state = stopped;
    private long elapsedTime = 0;
    private long startTime;

    private final JLabel time = new JLabel("time: 0");

    public void updateTime() {
        Date d = new Date();
        elapsedTime += (d.getTime() - startTime);
        time.setText(String.format("time: %d.%03d",
                                   elapsedTime/1000, elapsedTime%1000));
    }

    // Constructor
    //
    public UI(Coordinator c, MSTworld w,
            Animation a, long sd, RootPaneContainer pane) {
        final UI ui = this;
        coordinator = c;
        world = w;
        animation = a;

        final JPanel buttons = new JPanel();   // button panel
            final JButton runButton = new JButton("Run");
            final JButton pauseButton = new JButton("Pause");
            final JButton resetButton = new JButton("Reset");
            final JButton randomizeButton = new JButton("Randomize");
            final JButton quitButton = new JButton("Quit");

        final JPanel stats = new JPanel();   // statistics panel

        final JLabel seed = new JLabel("seed: " + sd + "   ");

        runButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                if (state == stopped) {
                    state = running;
                    root.setDefaultButton(pauseButton);
                    Worker w = new Worker(world, coordinator,
                                          ui, animation);
                    Date d = new Date();
                    startTime = d.getTime();
                    w.start();
                } else if (state == paused) {
                    state = running;
                    root.setDefaultButton(pauseButton);
                    Date d = new Date();
                    startTime = d.getTime();
                    coordinator.toggle();
                }
            }
        });
        pauseButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                if (state == running) {
                    updateTime();
                    state = paused;
                    root.setDefaultButton(runButton);
                    coordinator.toggle();
                }
            }
        });
        resetButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                state = stopped;
                coordinator.stop();
                root.setDefaultButton(runButton);
                world.reset();
                animation.reset();
                elapsedTime = 0;
                time.setText("time: 0");
            }
        });
        randomizeButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                state = stopped;
                coordinator.stop();
                root.setDefaultButton(runButton);
                long v = world.randomize();
                animation.reset();
                seed.setText("seed: " + v + "   ");
                elapsedTime = 0;
                time.setText("time: 0");
            }
        });
        quitButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                System.exit(0);
            }
        });

        // put the buttons into the button panel:
        buttons.setLayout(new FlowLayout());
        buttons.add(runButton);
        buttons.add(pauseButton);
        buttons.add(resetButton);
        buttons.add(randomizeButton);
        buttons.add(quitButton);

        // put the labels into the statistics panel:
        stats.add(seed);
        stats.add(time);

        // put the MSTworld canvas, the button panel, and the stats
        // label into the UI:
        setLayout(new BoxLayout(this, BoxLayout.Y_AXIS));
        setBorder(BorderFactory.createEmptyBorder(externalBorder,
            externalBorder, externalBorder, externalBorder));
        add(a);
        add(buttons);
        add(stats);

        // put the UI into the Frame
        pane.getContentPane().add(this);
        root = getRootPane();
        root.setDefaultButton(runButton);
    }
}
