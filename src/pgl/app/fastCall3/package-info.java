/**
 * The FastCall3 package provides high-performance tools for genetic variant calling and genotyping.
 * 
 * <p>This package implements the core functionality of the FastCall3 pipeline, which is designed for
 * efficient processing of high-throughput sequencing data to identify genetic variations and
 * determine genotypes across multiple samples.</p>
 *
 * <h2>Key Components:</h2>
 * <ul>
 *   <li>{@link FastCall3} - Main entry point for the FastCall3 pipeline</li>
 *   <li>{@link DiscoverVariationF3} - Discovers genetic variations from sequence data</li>
 *   <li>{@link BuildVariationLibraryF3} - Constructs a library of genetic variations</li>
 *   <li>{@link ViewVariationLibraryF3} - Views and converts variation libraries</li>
 *   <li>{@link ScanGenotypeF3} - Performs genotyping across samples</li>
 *   <li>{@link IndividualGenotypeF3} - Represents genotype data for a single individual</li>
 *   <li>{@link VariationLibraryF3} - Manages a collection of genetic variations</li>
 * </ul>
 *
 * <h2>Features:</h2>
 * <ul>
 *   <li>High-performance variant calling optimized for large datasets</li>
 *   <li>Memory-efficient data structures for handling large populations</li>
 *   <li>Support for both binary and text-based I/O formats</li>
 *   <li>Parallel processing capabilities for improved performance</li>
 *   <li>Flexible filtering and customization options</li>
 * </ul>
 *
 * <p>For more information, visit the
 * <a href="https://github.com/PlantGeneticsLab/TIGER">TIGER GitHub repository</a>.</p>
 *
 * @author Fei Lu
 * @version 3.0
 * @since 1.0
 */
package pgl.app.fastCall3;