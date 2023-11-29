
/**
 * This package provides the JMatrix matrix library. It intends to
 * <ul>
 *   <li>support garbage collection free operations;</li>
 *   <li>rely on minimal dependencies;</li>
 *   <li>and stay compilable on embedded Java implementations.</li>
 * </ul>
 *
 * Currently, only dense matrices are supported.<br/><br/>
 *
 * Supported matrix decompositions:
 * <ul>
 *   <li>LU decomposition;</li>
 *   <li>QR decomposition (full and reduced forms);</li>
 *   <li>SVD: singular value decomposition (full and reduced forms);</li>
 *   <li>Cholesky decomposition (for positive-definite matrices, LL and LDL forms).</li>
 * </ul>
 *
 * Other supported features:
 * <ul>
 *   <li>Constant-time matrix transpose;</li>
 *   <li>Extendable unary and binary elementwise operations;</li>
 *   <li>Orthogonalization by the Gram-Schmidt process;</li>
 *   <li>Solving equations by back substitution for triangular matrices;</li>
 *   <li>Manual placement of results in order to avoid garbage collection.</li>
 * </ul>
 */
package ai.gandg.jmatrix;
