#pragma once

namespace novac
{
/** Simple structure for storing the Id of a spectrometer.
    Spectrometers can normally not be identified by the serial only since there may be multiple
        spectrometer channels inside one device, hence is a channel number also needed. */
struct SpectrometerId
{
    SpectrometerId()
        : serial(""), channel(0)
    {
    }

    SpectrometerId(const std::string& serialNumber, int channelNumber)
        : serial(serialNumber), channel(channelNumber)
    {
    }

    /** The serial number of the device */
    std::string serial;

    /** The channel of the device, often zero but may be one for S2000 devices */
    int channel = 0;

    bool operator<(const SpectrometerId& other) const
    {
        int serialCompare = serial.compare(other.serial);
        if (serialCompare == 0)
        {
            return channel < other.channel;
        }
        else
        {
            return serialCompare < 0;
        }
    }

    bool operator==(const SpectrometerId& other) const
    {
        return serial.compare(other.serial) == 0 && channel == other.channel;
    }
};
}