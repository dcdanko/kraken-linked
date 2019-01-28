#include "seqreader.hpp"
#include "gtest/gtest.h"


namespace {


// Tests that the Foo::Bar() method does Abc.
TEST(SeqReaderTest, CanReadPlain) {
  const std::string input_filepath = "../tests/data/sample_fastq_zymo_controls";
  SeqReader *reader;
  reader = new FastqReader(input_filepath);
  ASSERT_EQ(
    reader.next_sequence().seq,
    'GGCATGACCTCTGCGAGGTCATGCCGCGGCCGGGCTGGGCGACGGGCTCGCCGATGCGGCGTCCCGCCCCCGCTGGTCCTCGGCGTCCCGGCCGGACGCGCTGGCGGAGCCCGAAGNGCTGGAGGAA'
  );
}

}  // namespace

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}