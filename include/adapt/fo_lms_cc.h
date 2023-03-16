/* -*- c++ -*- */
/*
 * Copyright 2023 gr-adapt author.
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */


#ifndef INCLUDED_ADAPT_FO_LMS_CC_H
#define INCLUDED_ADAPT_FO_LMS_CC_H

#include <adapt/api.h>
#include <gnuradio/block.h>

namespace gr {
namespace adapt {

/*!
 * \brief <+description of block+>
 * \ingroup adapt
 *
 */
class ADAPT_API fo_lms_cc : virtual public gr::block {
    protected:
    virtual gr_complex error(const gr_complex& desired, const gr_complex& out) = 0;
    virtual void update_tap(gr_complex& tap, const gr_complex& in) = 0;

    public:
    typedef std::shared_ptr<fo_lms_cc> sptr;

    /*!
     * \brief Return a shared_ptr to a new instance of adapt::fo_lms_cc.
     *
     * To avoid accidental use of raw pointers, adapt::fo_lms_cc's
     * constructor is in a private implementation
     * class. adapt::fo_lms_cc::make is the public interface for
     * creating new instances.
     */
    static sptr make(double samp_rate, int num_taps, float mu_cir, float mu_cfo, float mu_sfo);

    virtual float get_mu_cir() const = 0;
    virtual void set_mu_cir(float mu_cir) = 0;
    virtual float get_mu_cfo() const = 0;
    virtual void set_mu_cfo(float mu_cfo) = 0;
    virtual float get_mu_sfo() const = 0;
    virtual void set_mu_sfo(float mu_sfo) = 0;
};

} // namespace adapt
} // namespace gr

#endif /* INCLUDED_ADAPT_FO_LMS_CC_H */
