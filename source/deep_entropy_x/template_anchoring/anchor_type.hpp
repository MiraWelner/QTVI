#pragma once
//
// anchor_type.hpp
//
// AnchorType AND NOTHING ELSE. No includes, no dependencies, deliberately.
//
// The enum used to live in feature_marks.hpp, which pulls in template_bank.hpp
// (and through it noise_marking_gui/annotation_types.hpp), <functional> and the
// whole FeatureMarks class. Anything that merely wanted to NAME an alignment --
// landmark_admissibility.hpp, anchor_view.hpp, template_io.hpp -- had to drag
// all of that in with it, which is why:
//
//   * landmark_admissibility.hpp carries a warning that feature_marks.hpp must
//     not include it back "or the graph cycles", and asks in as many words for
//     this header to exist: "If AnchorType is ever split into its own
//     dependency-free header, switch this to that instead."
//   * template_io.hpp types its raw_anchors map as `int` rather than
//     AnchorType, with a note explaining that it settled for int because
//     including feature_marks.hpp there was not acceptable. It can now be
//     typed properly, though this change does not do that -- it is a format-
//     adjacent edit and belongs on its own.
//
// feature_marks.hpp includes this and re-exports the name, so every existing
// user of AnchorType keeps compiling unchanged.
//
// J_POINT == S_END: the alignment anchored on the J point. anchor_view.hpp
// labels it "T" for the operator, because making the ST-T segment sharp is
// what that alignment is for -- see the note there.
//
enum class AnchorType { P_ONSET, Q_ONSET, R_PEAK, J_POINT };