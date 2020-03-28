package pgl.infra.dna;

/**
 * FastaRecord with 3 bits storage
 * @author feilu
 */
public class FastaRecordBit extends Sequence3Bit implements FastaRecordInterface {
    String name;
    int id;

    public FastaRecordBit(String name, String seq, int id) {
        super(seq);
        this.name = name;
        this.id = id;
    }
    
    public FastaRecordBit (String name, Sequence3Bit s3, int id) {
        super.sequenceLength = s3.sequenceLength;
        super.seqS = s3.seqS;
        this.name = name;
        this.id = id;
    }
    
    @Override
    public String getName() {
        return name;
    }

    @Override
    public int getID() {
        return id;
    }

    @Override
    public void setName(String newName) {
        name = newName;
    }

    @Override
    public void setID(int id) {
        this.id = id;
    }
}
