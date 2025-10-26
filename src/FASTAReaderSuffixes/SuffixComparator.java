package FASTAReaderSuffixes;

import java.util.Comparator;

/**
 * A custom Comparator for sorting Suffix objects based on the content of the
 * FASTA file. It needs access to the FASTA file's data array.
 *
 * @author mmiguel, rgarciacarmona
 */
public class SuffixComparator implements Comparator<Suffix> {

	// Since FASTAReader is the base class containing content and validBytes,
	// we need a reference to it to access the data.
	private final FASTAReaderSuffixes reader;

	/**
	 * Constructs a SuffixComparator, initializing it with the data provider.
	 * 
	 * @param reader The FASTAReader object that holds the FASTA sequence data.
	 */
	public SuffixComparator(FASTAReaderSuffixes reader) {
		this.reader = reader;
	}

	/**
	 * Compares two Suffix objects by comparing the byte sequences starting at their
	 * respective indices in the FASTA data array.
	 * 
	 * @param s1 The first Suffix.
	 * @param s2 The second Suffix.
	 * @return -1, 0, or 1 as the first suffix is less than, equal to, or greater
	 *         than the second.
	 */
	@Override
	public int compare(Suffix s1, Suffix s2) {
		int index1 = s1.suffixIndex;
		int index2 = s2.suffixIndex;

		int len = Math.min(this.reader.getValidBytes() - index1, this.reader.getValidBytes() - index2);

		for (int i = 0; i < len; i++) {
			if (this.reader.getContent()[index1 + i] < this.reader.getContent()[index2 + i])
				return -1;
			if (this.reader.getContent()[index1 + i] > this.reader.getContent()[index2 + i])
				return 1;
		}
		// If one is a prefix of the other, the shorter one (with larger index) comes
		// first (lexicographically).
		return index2 - index1;
	}
}