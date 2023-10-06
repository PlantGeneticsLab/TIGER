package pgl.infra.dna;

/**
 * Fasta record with byte storage.
 * @author feilu
 */
public class FastaRecordByte extends SequenceByte implements FastaRecordInterface {
    /**
     * Name of the Fasta record.
     */
    String description;
    /**
     * ID of the Fasta record.
     */
    int id;

    /**
     * Constructs an object.
     * @param description
     * @param seq
     * @param id
     */
    public FastaRecordByte(String description, String seq, int id) {
        super(seq);
        this.description = description;
        this.id = id;
    }

    /**
     * Constructs an object with {@link SequenceByte}.
     * @param description
     * @param sb
     * @param id
     */
    public FastaRecordByte(String description, SequenceByte sb, int id) {
        super.seqAscII = sb.seqAscII;
        this.description = description;
        this.id = id;
    }
    
    @Override
    public String getDescription() {
        return description;
    }

    @Override
    public int getID() {
        return id;
    }

    @Override
    public void setDescription(String description) {
        this.description = description;
    }

    @Override
    public void setID(int id) {
        this.id = id;
    }

}
