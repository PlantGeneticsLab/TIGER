package pgl.infra.dna;

/**
 * FastaRecord with byte storage
 * @author feilu
 */
public class FastaRecordByte extends SequenceByte implements FastaRecordInterface {
    String name;
    int id;

    public FastaRecordByte(String name, String seq, int id) {
        super(seq);
        this.name = name;
        this.id = id;
    }
    
    public FastaRecordByte(String name, SequenceByte sb, int id) {
        super.seqAscII = sb.seqAscII;
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
