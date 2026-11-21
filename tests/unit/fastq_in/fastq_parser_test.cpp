// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2026 Mikkel Schubert <mikkelsch@gmail.com>
#include "../scheduler_test.hpp" // for STEP_PROCESS_TO, ...
#include "commontypes.hpp"       // for output_format, read_mate, read_type
#include "errors.hpp"            // for assert_failed, fastq_error
#include "fastq_in.hpp"          // for fastq_parser, reader_chunk, ...
#include "fastq_in_test.hpp"     // for DEFAULT_BUFFER_CAPACITY, ...
#include "logging.hpp"           // for log_capture
#include "scheduler.hpp"         // buffer, chunk_ptr, fastq, ...
#include "serializer.hpp"        // for read_meta, serializer
#include "testing.hpp"           // for GENERATE, INFO, REQUIRE, ...
#include "testutils.hpp"         // for random_fastq_buffer
#include "utilities.hpp"         // for dynamic_cast_unique
#include <limits>                // for numeric_limits
#include <memory>                // for shared_ptr
#include <tuple>                 // for ignore

namespace adapterremoval {

#define PARAMETERIZE_INPUT_TYPE                                                \
  GENERATE(fastq_input_type::mate_1,                                           \
           fastq_input_type::mate_2,                                           \
           fastq_input_type::interleaved)

#define PARAMETERIZE_INPUT_12                                                  \
  GENERATE(fastq_input_type::mate_1, fastq_input_type::mate_2)

namespace {

constexpr size_t DEFAULT_HEAD = std::numeric_limits<size_t>::max();

class parser_builder
{
public:
  parser_builder()
    : m_state(make_state())
  {
  }

  parser_builder(std::shared_ptr<reader_state> state)
    : m_state(std::move(state))
  {
  }

  parser_builder& set_input_type(fastq_input_type value)
  {
    m_input_type = value;
    return *this;
  }

  parser_builder& set_next_step(step_id value)
  {
    m_next_step = value;
    return *this;
  }

  parser_builder& set_head(size_t value)
  {
    m_head = value;
    return *this;
  }

  parser_builder& set_filenames(std::vector<std::string> value)
  {
    m_filenames = std::move(value);
    return *this;
  }

  [[nodiscard]] fastq_parser build() const
  {
    return fastq_parser{ m_filenames,
                         m_input_type,
                         m_head,
                         m_state,
                         m_next_step };
  }

private:
  std::shared_ptr<reader_state> m_state;
  fastq_input_type m_input_type = fastq_input_type::mate_1;
  step_id m_next_step = DEFAULT_NEXT_STEP;
  size_t m_head = DEFAULT_HEAD;
  std::vector<std::string> m_filenames{ "filename1" };
};

reader_chunk_ptr
make_uncompressed_chunk(std::shared_ptr<reader_state>& state,
                        chunk_eof eof,
                        std::string_view data = {})
{
  auto chunk = make_reader_chunk(state, eof, reader_chunk::type::text);
  chunk->data.append(data);
  return chunk;
}

void
fastqs_into_buffer(const std::vector<fastq>& fastqs, buffer& dst)
{
  read_meta m{ read_type::se };
  serializer s{ output_format::fastq };
  for (const auto& it : fastqs) {
    s.record(dst, it, m);
  }
}

std::vector<fastq>
gather_reads(const fastq_parser_chunk& chunk, fastq_input_type type)
{
  std::vector<fastq> output;
  switch (type) {
    case fastq_input_type::mate_1:
      REQUIRE((chunk.type == fastq_input_type::mate_1));
      // NOLINTNEXTLINE(readability-container-size-empty)
      REQUIRE(chunk.reads_2 == std::vector<fastq>{});
      output = chunk.reads_1;
      break;
    case fastq_input_type::mate_2:
      REQUIRE((chunk.type == fastq_input_type::mate_2));
      // NOLINTNEXTLINE(readability-container-size-empty)
      REQUIRE(chunk.reads_1 == std::vector<fastq>{});
      output = chunk.reads_2;
      break;
    case fastq_input_type::interleaved: {
      REQUIRE((chunk.type == fastq_input_type::interleaved));
      REQUIRE(chunk.reads_1.size() == chunk.reads_2.size());
      for (size_t i = 0; i < chunk.reads_1.size(); ++i) {
        output.emplace_back(chunk.reads_1.at(i));
        output.emplace_back(chunk.reads_2.at(i));
      }
      break;
    }
    default:
      AR_FAIL("invalid input_type");
  }

  return output;
}

} // namespace

////////////////////////////////////////////////////////////////////////////////

TEST_CASE("fastq_parser initialization", "[fastq_parser]")
{
  auto step = parser_builder().build();
  REQUIRE(step.next_step() == DEFAULT_NEXT_STEP);
  REQUIRE(step.name() == "fastq_parser");
  REQUIRE(step.ordering() == processing_order::ordered);
}

TEST_CASE("fastq_parser throws on empty filenames list", "[fastq_parser]")
{
  auto builder = parser_builder() //
                   .set_filenames({})
                   .set_input_type(PARAMETERIZE_INPUT_TYPE);

  REQUIRE_THROWS_AS(builder.build(), assert_failed);
}

TEST_CASE("fastq_parser throws on nullptr state", "[fastq_parser]")
{
  auto builder = parser_builder(nullptr) //
                   .set_input_type(PARAMETERIZE_INPUT_TYPE);

  REQUIRE_THROWS_CONTAINS(builder.build(), assert_failed, "~5a42dfdc29f8abdb~");
}

TEST_CASE("fastq_parser returns empty vector for file with no non-empty lines",
          "[fastq_parser]")
{
  log::log_capture cap;
  auto state = make_state();
  auto step = parser_builder{ state } //
                .set_input_type(PARAMETERIZE_INPUT_TYPE)
                .build();

  auto input = make_uncompressed_chunk(state, chunk_eof::yes);
  input->data.append(GENERATE("", "\n\n\n", "\r\n\r\n\r\n", "\r\n\n\r\n"));

  std::ignore = STEP_SIGNAL_TO(step_id, step);
  auto output =
    STEP_PROCESS_TO_CHUNK(fastq_parser_chunk, step, std::move(input));
  // NOLINTNEXTLINE(readability-container-size-empty)
  REQUIRE(output->reads_1 == std::vector<fastq>{});
  // NOLINTNEXTLINE(readability-container-size-empty)
  REQUIRE(output->reads_2 == std::vector<fastq>{});
  REQUIRE(output->eof);

  REQUIRE_NOTHROW(state->buffers.finalize());
}

TEST_CASE("fastq_parser returns EOF only when getline() returns EOF",
          "[fastq_parser]")
{
  const auto input_type = PARAMETERIZE_INPUT_TYPE;
  const auto fastq_size = (input_type == fastq_input_type::interleaved) ? 1 : 2;
  auto state = make_state(DEFAULT_BUFFER_CAPACITY, fastq_size);
  auto step = parser_builder(state) //
                .set_input_type(input_type)
                .build();

  auto input = make_uncompressed_chunk(state,
                                       chunk_eof::yes,
                                       "@read1\n"
                                       "GTGAATTAGGAATGTCAGCA\n"
                                       "+\n"
                                       "=;>9@A?8A:<;@?<==7?>\n"
                                       "@read2\n"
                                       "TAGTTATTTTGACGTACTGG\n"
                                       "+\n"
                                       "98456;58;=8<::4=58<7");
  input->data.append(GENERATE("", "\n", "\r\n", "\n\n"));

  std::ignore = STEP_SIGNAL_TO(step_id, step);
  auto output =
    STEP_PROCESS_TO_CHUNK(fastq_parser_chunk, step, std::move(input));
  REQUIRE(gather_reads(*output, input_type) ==
          std::vector<fastq>{
            fastq("@read1", "GTGAATTAGGAATGTCAGCA", "=;>9@A?8A:<;@?<==7?>"),
            fastq("@read2", "TAGTTATTTTGACGTACTGG", "98456;58;=8<::4=58<7"),
          });
  REQUIRE_FALSE(output->eof);

  output = STEP_SIGNAL_TO_CHUNK(fastq_parser_chunk, step);
  // NOLINTNEXTLINE(readability-container-size-empty)
  REQUIRE(output->reads_1 == std::vector<fastq>{});
  // NOLINTNEXTLINE(readability-container-size-empty)
  REQUIRE(output->reads_2 == std::vector<fastq>{});
  REQUIRE(output->eof);

  log::log_capture cap;
  REQUIRE_NOTHROW(state->buffers.finalize());
}

TEST_CASE("fastq_parser validates basic FASTQ properties", "[fastq_parser]")
{
  const auto input_type = PARAMETERIZE_INPUT_12;
  auto state = make_state();
  auto step = parser_builder(state) //
                .set_input_type(input_type)
                .build();
  auto input = make_uncompressed_chunk(state, chunk_eof::yes);

  std::ignore = STEP_SIGNAL_TO(step_id, step);

  SECTION("proper FASTQ header required")
  {
    input->data.append(GENERATE(">read1\n", "@\n", "\n"));
    input->data.append("GTGAATTAGGAATGTCAGCA\n"
                       "+\n"
                       "=;>9@A?8A:<;@?<==7?>\n");

    REQUIRE_THROWS_CONTAINS(step.process(std::move(input)),
                            fastq_error,
                            "Malformed or empty FASTQ header");
  }

  SECTION("sequence must be non-empty")
  {
    input->data.append("@read1\n"
                       "\n"
                       "+\n"
                       "=;>9@A?8A:<;@?<==7?>\n");

    REQUIRE_THROWS_CONTAINS(step.process(std::move(input)),
                            fastq_error,
                            "sequence is empty");
  }

  SECTION("separator must be present")
  {
    input->data.append("@read1\n"
                       "GTGAATTAGGAATGTCAGCA");
    input->data.append(GENERATE("-\n", "\n"));
    input->data.append("=;>9@A?8A:<;@?<==7?>\n");

    REQUIRE_THROWS_CONTAINS(step.process(std::move(input)),
                            fastq_error,
                            "FASTQ record lacks separator character");
  }

  SECTION("text after separator is ignored")
  {
    input->data.append("@read1\n"
                       "GTGAATTAGGAATGTCAGCA\n");
    input->data.append(GENERATE("+\n", "+read1\n", "+random text\n"));
    input->data.append("=;>9@A?8A:<;@?<==7?>\n");

    auto chunk =
      STEP_PROCESS_TO_CHUNK(fastq_parser_chunk, step, std::move(input));
    REQUIRE(gather_reads(*chunk, input_type) ==
            std::vector<fastq>{
              { "@read1", "GTGAATTAGGAATGTCAGCA", "=;>9@A?8A:<;@?<==7?>" } });
  }

  SECTION("quality length must match sequence length")
  {
    input->data.append("@read1\n"
                       "GTGAATTAGGAATGTCAGCA\n"
                       "+\n");
    input->data.append(
      GENERATE("=;>9@A?8A:<;@?<==7?\n", "=;>9@A?8A:<;@?<==7?>!\n", "\n"));

    REQUIRE_THROWS_CONTAINS(step.process(std::move(input)),
                            fastq_error,
                            "sequence/quality lengths do not match");
  }

  SECTION("sequence/quality content is not validated")
  {
    input->data.append("@read1\n"
                       "This isn't real DNA!\n"
                       "+\n"
                       "=;>9\tA?8A:<;@?<= 7?\x7f\n");

    auto chunk =
      STEP_PROCESS_TO_CHUNK(fastq_parser_chunk, step, std::move(input));
    auto reads = gather_reads(*chunk, input_type);
    REQUIRE(reads.size() == 1);
    REQUIRE(reads.front().header() == "@read1");
    REQUIRE(reads.front().sequence() == "This isn't real DNA!");
    REQUIRE(reads.front().qualities() == "=;>9\tA?8A:<;@?<= 7?\x7f");
  }
}

TEST_CASE("fastq_parser throws if records are truncated", "[fastq_parser]")
{
  auto state = make_state();
  auto step = parser_builder(state) //
                .set_input_type(PARAMETERIZE_INPUT_TYPE)
                .build();
  auto input = make_uncompressed_chunk(state, chunk_eof::yes);

  if (GENERATE(true, false)) {
    input->data.append("@read1\n"
                       "GTGAATTAGGAATGTCAGCA\n"
                       "+\n"
                       "=;>9@A?8A:<;@?<==7?>\n");
  }

  std::ignore = STEP_SIGNAL_TO(step_id, step);

  SECTION("truncated after header")
  {
    input->data.append("@read2\n");
    REQUIRE_THROWS_CONTAINS(step.process(std::move(input)),
                            fastq_error,
                            "partial FASTQ record; cut off after header");
  }

  SECTION("truncated after sequence")
  {
    input->data.append("@read2\n"
                       "TAGTTATTTTGACGTACTGG\n");
    REQUIRE_THROWS_CONTAINS(step.process(std::move(input)),
                            fastq_error,
                            "partial FASTQ record; cut off after sequence");
  }

  SECTION("truncated after separator")
  {
    input->data.append("@read2\n"
                       "TAGTTATTTTGACGTACTGG\n"
                       "+\n");
    REQUIRE_THROWS_CONTAINS(step.process(std::move(input)),
                            fastq_error,
                            "partial FASTQ record; cut off after separator");
  }
}

TEST_CASE("fastq_parser handles both CR, LFCR, and mixes", "[fastq_parser]")
{
  const auto input_type = PARAMETERIZE_INPUT_12;
  auto state = make_state();
  auto step = parser_builder(state) //
                .set_input_type(input_type)
                .build();

  auto input = make_uncompressed_chunk(state, chunk_eof::yes);
  input->data.append("@read1");
  input->data.append(GENERATE("\n", "\r\n"));
  input->data.append("GTGAATTAGGAATGTCAGCA");
  input->data.append(GENERATE("\n", "\r\n"));
  input->data.append("+");
  input->data.append(GENERATE("\n", "\r\n"));
  input->data.append("=;>9@A?8A:<;@?<==7?>");
  input->data.append(GENERATE("", "\n", "\r\n"));

  std::ignore = STEP_SIGNAL_TO(step_id, step);
  auto output =
    STEP_PROCESS_TO_CHUNK(fastq_parser_chunk, step, std::move(input));
  REQUIRE(gather_reads(*output, input_type) ==
          std::vector<fastq>{
            { "@read1", "GTGAATTAGGAATGTCAGCA", "=;>9@A?8A:<;@?<==7?>" } });
  REQUIRE(output->eof);
}

TEST_CASE("fastq_parser returns nothing/suspend if parsing is done",
          "[fastq_parser]")
{
  auto state = make_state();
  auto step = parser_builder(state).build();
  state->mark_parsing(GENERATE(reader_state::parsing::reached_head,
                               reader_state::parsing::reached_eof));

  SECTION("signal")
  {
    auto signal = STEP_SIGNAL_TO(scheduler_signal, step);
    REQUIRE(signal == scheduler_signal::suspend);
  }

  SECTION("chunk")
  {
    auto input = make_uncompressed_chunk(state,
                                         chunk_eof::yes,
                                         "@read1\n"
                                         "GTGAATTAGGAATGTCAGCA\n"
                                         "+\n"
                                         "=;>9@A?8A:<;@?<==7?>\n");
    // fastq_parser may receive "dropped" chunks if state->done
    input->data_type = GENERATE(reader_chunk::type::unknown,
                                reader_chunk::type::bgzf,
                                reader_chunk::type::text);

    auto chunks = step.process(std::move(input));
    // NOLINTNEXTLINE(readability-container-size-empty)
    REQUIRE(chunks == chunk_vec{});
  }
}

TEST_CASE("fastq_parser requires non-compressed data", "[fastq_parser]")
{
  auto state = make_state();
  auto step = parser_builder(state).build();

  auto input = make_uncompressed_chunk(state,
                                       chunk_eof::yes,
                                       "@read1\n"
                                       "GTGAATTAGGAATGTCAGCA\n"
                                       "+\n"
                                       "=;>9@A?8A:<;@?<==7?>\n");

  SECTION("non-compressed")
  {
    input->data_type = reader_chunk::type::text;
    REQUIRE_NOTHROW(step.process(std::move(input)));
  }

  SECTION("unknown or compressed")
  {
    input->data_type =
      GENERATE(reader_chunk::type::unknown, reader_chunk::type::bgzf);
    REQUIRE_THROWS_CONTAINS(step.process(std::move(input)),
                            assert_failed,
                            "~a12ddd10369e1b5c~");
  }
}

TEST_CASE("fastq_parser requires non-nullptr chunk", "[fastq_parser]")
{
  auto state = make_state();
  auto step = parser_builder(state).build();
  if (GENERATE(true, false)) {
    state->mark_parsing(GENERATE(reader_state::parsing::reached_head,
                                 reader_state::parsing::reached_eof));
  }

  REQUIRE_THROWS_CONTAINS(step.process(nullptr),
                          assert_failed,
                          "~8ed559b0c3086856~");
}

TEMPLATE_TEST_CASE("fastq_parser accepts only reader_chunks",
                   "[fastq_parser]",
                   adapterremoval::fastq_parser_chunk,
                   adapterremoval::output_chunk)
{
  auto state = make_state();
  auto step = parser_builder(state).build();
  if (GENERATE(true, false)) {
    state->mark_parsing(GENERATE(reader_state::parsing::reached_head,
                                 reader_state::parsing::reached_eof));
  }

  auto chunk = std::make_unique<TestType>();
  REQUIRE_THROWS_CONTAINS(step.process(std::move(chunk)),
                          assert_failed,
                          "~8ed559b0c3086856~");
}

TEST_CASE("fastq_parser does not require newline between files",
          "[fastq_parser]")
{
  const auto input_type = PARAMETERIZE_INPUT_12;
  auto state = make_state(DEFAULT_BUFFER_CAPACITY, 1);
  auto step = parser_builder(state) //
                .set_filenames({ "filename1", "filename2" })
                .set_input_type(input_type)
                .build();

  auto input = make_uncompressed_chunk(state,
                                       chunk_eof::yes,
                                       "@read1\n"
                                       "GTGAATTAGGAATGTCAGCA\n"
                                       "+\n"
                                       "=;>9@A?8A:<;@?<==7?>\n"
                                       "@read2\n"
                                       "TAGTTATTTTGACGTACTGG\n"
                                       "+\n"
                                       "98456;58;=8<::4=58<7");
  input->data.append(GENERATE("", "\n", "\n\n"));

  std::ignore = STEP_SIGNAL_TO(step_id, step);
  auto output =
    STEP_PROCESS_TO_CHUNK(fastq_parser_chunk, step, std::move(input));
  REQUIRE(gather_reads(*output, input_type) ==
          std::vector<fastq>{
            { "@read1", "GTGAATTAGGAATGTCAGCA", "=;>9@A?8A:<;@?<==7?>" } });

  input = make_uncompressed_chunk(state,
                                  chunk_eof::yes,
                                  "@read3\n"
                                  "GTGTCGAGTTTGCTTATTGC\n"
                                  "+\n"
                                  "?;7@<;<??7@=:=>A<==7\n");

  output = STEP_SIGNAL_TO_CHUNK(fastq_parser_chunk, step);
  REQUIRE(gather_reads(*output, input_type) ==
          std::vector<fastq>{
            { "@read2", "TAGTTATTTTGACGTACTGG", "98456;58;=8<::4=58<7" } });

  std::ignore = STEP_SIGNAL_TO(step_id, step);
  output = STEP_PROCESS_TO_CHUNK(fastq_parser_chunk, step, std::move(input));
  REQUIRE(gather_reads(*output, input_type) ==
          std::vector<fastq>{
            { "@read3", "GTGTCGAGTTTGCTTATTGC", "?;7@<;<??7@=:=>A<==7" } });
}

TEST_CASE("fastq_parser::process returns chunk if able, signal otherwise",
          "[fastq_parser]")
{
  const auto input_type = PARAMETERIZE_INPUT_12;
  auto state = make_state(DEFAULT_BUFFER_CAPACITY, 32);
  auto step = parser_builder(state) //
                .set_input_type(input_type)
                .build();
  buffer total_input;
  buffer total_output;

  const auto signal_target = make_step_id_for_testing(7913);
  set_signal_target_for_testing(step, signal_target, read_mate::_1);

  std::ignore = STEP_SIGNAL_TO(step_id, step);
  for (const auto size : std::array{ 15, 16 }) {
    auto input = make_uncompressed_chunk(state, chunk_eof::no);
    input->data.append(random_fastq_buffer(size));
    total_input.append(input->data);

    auto signal = STEP_PROCESS_TO(step_id, step, std::move(input));
    REQUIRE(signal == signal_target);
    signal = STEP_SIGNAL_TO(step_id, step);
    REQUIRE(signal == signal_target);
  }

  auto input = make_uncompressed_chunk(state, chunk_eof::no);
  input->data.append(random_fastq_buffer(7));
  total_input.append(input->data);
  auto chunk =
    STEP_PROCESS_TO_CHUNK(fastq_parser_chunk, step, std::move(input));
  auto reads = gather_reads(*chunk, input_type);
  fastqs_into_buffer(reads, total_output);
  REQUIRE(reads.size() == 32);

  input = make_uncompressed_chunk(state, chunk_eof::yes);
  chunk = STEP_PROCESS_TO_CHUNK(fastq_parser_chunk, step, std::move(input));
  reads = gather_reads(*chunk, input_type);
  REQUIRE(reads.size() == 6);

  fastqs_into_buffer(reads, total_output);
  REQUIRE(total_input == total_output);

  auto ssignal = STEP_SIGNAL_TO(scheduler_signal, step);
  REQUIRE(ssignal == scheduler_signal::suspend);

  log::log_capture cap;
  REQUIRE_NOTHROW(state->buffers.finalize());
}

TEST_CASE("fastq_parser::signal returns chunk if able, signal otherwise",
          "[fastq_parser]")
{
  const auto input_type = PARAMETERIZE_INPUT_12;
  auto state = make_state(DEFAULT_BUFFER_CAPACITY, 16);
  auto step = parser_builder(state) //
                .set_input_type(input_type)
                .build();
  buffer total_input;
  buffer total_output;

  const auto signal_target = make_step_id_for_testing(7913);
  set_signal_target_for_testing(step, signal_target, read_mate::_1);

  auto input = make_uncompressed_chunk(state, chunk_eof::no);
  input->data.append(random_fastq_buffer(49));
  total_input.append(input->data);

  STEP_PROCESS_TO_NOTHING(step, std::move(input));
  auto chunk = STEP_SIGNAL_TO_CHUNK(fastq_parser_chunk, step);
  auto reads = gather_reads(*chunk, input_type);
  fastqs_into_buffer(reads, total_output);
  REQUIRE(reads.size() == 16);

  chunk = STEP_SIGNAL_TO_CHUNK(fastq_parser_chunk, step);
  reads = gather_reads(*chunk, input_type);
  fastqs_into_buffer(reads, total_output);
  REQUIRE(reads.size() == 16);
  chunk = STEP_SIGNAL_TO_CHUNK(fastq_parser_chunk, step);
  reads = gather_reads(*chunk, input_type);
  fastqs_into_buffer(reads, total_output);
  REQUIRE(reads.size() == 16);

  auto signal = STEP_SIGNAL_TO(step_id, step);
  REQUIRE(signal == signal_target);

  input = make_uncompressed_chunk(state, chunk_eof::yes);
  chunk = STEP_PROCESS_TO_CHUNK(fastq_parser_chunk, step, std::move(input));
  reads = gather_reads(*chunk, input_type);
  fastqs_into_buffer(reads, total_output);
  REQUIRE(reads.size() == 1);

  auto ssignal = STEP_SIGNAL_TO(scheduler_signal, step);
  REQUIRE(ssignal == scheduler_signal::suspend);

  log::log_capture cap;
  REQUIRE(total_input == total_output);
  REQUIRE_NOTHROW(state->buffers.finalize());
}

TEST_CASE("fastq_parser returns multiple pending chunks", "[fastq_parser]")
{
  auto state = make_state(DEFAULT_BUFFER_CAPACITY, 100);
  auto step = parser_builder(state) //
                .build();

  const size_t pending = GENERATE(range(1, 5));
  for (size_t i = 0; i < pending; ++i) {
    std::ignore = STEP_SIGNAL_TO(step_id, step);
  }

  auto fastqs = random_fastq_buffer(399);
  auto input_chunk = make_uncompressed_chunk(state, chunk_eof::yes);
  input_chunk->data.append(fastqs);

  buffer total_output;
  auto chunks = step.process(std::move(input_chunk));
  REQUIRE(chunks.size() == pending);
  for (auto& chunk : chunks) {
    auto* ptr = std::get_if<chunk_ptr>(&chunk);
    REQUIRE(ptr);
    auto output = dynamic_cast_unique<fastq_parser_chunk>(*ptr);
    REQUIRE(output.get());

    auto reads = gather_reads(*output, fastq_input_type::mate_1);
    REQUIRE((reads.size() == 100 || (output->eof && reads.size() == 99)));
    fastqs_into_buffer(reads, total_output);
  }

  REQUIRE(total_output.size() <= fastqs.size());
  buffer prefix;
  prefix.append(fastqs.data(), total_output.size());
  REQUIRE(total_output == prefix);
}

TEST_CASE("fastq_parser::finalize requires EOF or reached head",
          "[fastq_parser]")
{
  auto state = make_state();
  auto step = parser_builder(state).build();

  SECTION("neither at head nor EOF")
  {
    REQUIRE_THROWS_CONTAINS(step.finalize(),
                            assert_failed,
                            "~1939d44ea0df7a34~");
  }

  SECTION("only parser reached eof (impossible state)")
  {
    state->mark_parsing(reader_state::parsing::reached_eof);
    REQUIRE_THROWS_CONTAINS(step.finalize(),
                            assert_failed,
                            "~1939d44ea0df7a34~");
  }

  SECTION("parser reached head")
  {
    state->mark_parsing(reader_state::parsing::reached_head);
    REQUIRE_NOTHROW(step.finalize());
  }

  SECTION("parser reached EOF")
  {
    auto input = make_uncompressed_chunk(state, chunk_eof::yes);
    std::ignore = STEP_SIGNAL_TO(step_id, step);
    std::ignore =
      STEP_PROCESS_TO_CHUNK(fastq_parser_chunk, step, std::move(input));

    REQUIRE_FALSE(state->is_parsing());
    REQUIRE_NOTHROW(step.finalize());
    REQUIRE_NOTHROW(state->buffers.finalize());
  }
}

TEST_CASE("fastq_parser::finalize releases buffers if parsing is done",
          "[fastq_parser]")
{
  auto state = make_state();
  auto step = parser_builder(state).build();

  auto input = make_uncompressed_chunk(state, chunk_eof::yes);
  input->data.append(random_fastq_buffer(15));
  STEP_PROCESS_TO_NOTHING(step, std::move(input));

  log::log_capture cap;
  REQUIRE_THROWS_CONTAINS(state->buffers.finalize(),
                          assert_failed,
                          "not all objects returned to buffers cache");

  state->mark_parsing(reader_state::parsing::reached_head);
  step.finalize();
  REQUIRE_NOTHROW(state->buffers.finalize());
}

TEST_CASE("fastq_parser joins chunks from the same file", "[fastq_parser]")
{
  const auto input_type = PARAMETERIZE_INPUT_TYPE;
  const std::string_view input_data = "@read1\n"
                                      "GTGAATTAGGAATGTCAGCA\n"
                                      "+\n"
                                      "=;>9@A?8A:<;@?<==7?>\n"
                                      "@read2\n"
                                      "TAGTTATTTTGACGTACTGG\n"
                                      "+\n"
                                      "98456;58;=8<::4=58<7\n";

  for (size_t i = 0; i <= input_data.size(); i++) {
    INFO("i = " << i << " of " << input_data.size());
    auto state = make_state();
    auto step = parser_builder(state) //
                  .set_input_type(input_type)
                  .build();
    buffer output;

    auto input = make_uncompressed_chunk(state, chunk_eof::no);
    input->data.append(input_data.substr(0, i));
    std::ignore = STEP_SIGNAL_TO(step_id, step);
    std::ignore = STEP_PROCESS_TO(step_id, step, std::move(input));

    input = make_uncompressed_chunk(state, chunk_eof::yes);
    input->data.append(input_data.substr(i));
    std::ignore = STEP_SIGNAL_TO(step_id, step);
    auto chunk =
      STEP_PROCESS_TO_CHUNK(fastq_parser_chunk, step, std::move(input));
    REQUIRE(chunk->eof);

    fastqs_into_buffer(gather_reads(*chunk, input_type), output);

    buffer expected_output;
    expected_output.append(input_data.data(), input_data.size());
    REQUIRE(output == expected_output);

    log::log_capture cap;
    REQUIRE_NOTHROW(state->buffers.finalize());
  }
}

TEST_CASE("fastq_parser respects head", "[fastq_parser]")
{
  const size_t reads = 32;
  const size_t head = GENERATE(0, 1, 2, 31, 32);
  const auto input_type = PARAMETERIZE_INPUT_TYPE;
  const auto multiplier = input_type == fastq_input_type::interleaved ? 2 : 1;
  std::vector<std::string> filenames{ "filename1" };
  if (GENERATE(true, false)) {
    filenames.emplace_back("filename2");
  }

  auto state = make_state();
  auto step = parser_builder(state) //
                .set_input_type(input_type)
                .set_filenames(filenames)
                .set_head(head)
                .build();

  auto chunk =
    make_uncompressed_chunk(state, GENERATE(chunk_eof::yes, chunk_eof::no));
  chunk->data.append(random_fastq_buffer(reads * multiplier));
  STEP_PROCESS_TO_NOTHING(step, std::move(chunk));
  auto output = STEP_SIGNAL_TO_CHUNK(fastq_parser_chunk, step);
  REQUIRE(gather_reads(*output, input_type).size() == head * multiplier);
  REQUIRE(output->eof);
}

TEST_CASE("fastq_parser respects head past EOF", "[fastq_parser]")
{
  const size_t reads = 32;
  const size_t head = GENERATE(32, 33, 64);
  const auto input_type = PARAMETERIZE_INPUT_TYPE;
  const auto multiplier = input_type == fastq_input_type::interleaved ? 2 : 1;
  std::vector<std::string> filenames{ "filename1" };

  auto state = make_state();
  auto step = parser_builder(state) //
                .set_input_type(input_type)
                .set_filenames(filenames)
                .set_head(head)
                .build();

  auto chunk = make_uncompressed_chunk(state, chunk_eof::yes);
  chunk->data.append(random_fastq_buffer(reads * multiplier));
  std::ignore = STEP_SIGNAL_TO(step_id, step);
  auto output =
    STEP_PROCESS_TO_CHUNK(fastq_parser_chunk, step, std::move(chunk));
  REQUIRE(gather_reads(*output, input_type).size() == reads * multiplier);
  REQUIRE(output->eof);
}

TEST_CASE("fastq_parser catches uneven interleaved reads at head",
          "[fastq_parser]")
{
  const size_t total_pairs = GENERATE(1, 2, 3);
  auto state = make_state();
  auto step = parser_builder(state) //
                .set_input_type(fastq_input_type::interleaved)
                .set_head(total_pairs)
                .build();

  auto chunk = make_uncompressed_chunk(state, chunk_eof::yes);
  chunk->data.append(random_fastq_buffer((total_pairs * 2) - 1));

  std::ignore = STEP_SIGNAL_TO(step_id, step);
  REQUIRE_THROWS_CONTAINS(step.process(std::move(chunk)),
                          fastq_error,
                          "Found uneven number of interleaved reads");
}

} // namespace adapterremoval
