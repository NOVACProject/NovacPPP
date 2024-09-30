#pragma once

#include <PPPLib/MFC/CString.h>
#include <SpectralEvaluation/Log.h>

/** Appends an information / warning message to the logs.
    Implementations are in Common.cpp (using Poco) */
void ShowMessage(const novac::CString& message);
void ShowMessage(const char message[]);
void ShowMessage(const novac::CString& message, novac::CString connectionID);
void ShowMessage(const std::string& message);

/** Appends an error message to the logs */
void ShowError(const novac::CString& message);
void ShowError(const char message[]);
