/*
 *
 * MetaTarget NamedThreadFactory
 *
 * Copyright (C) 2022 Anestis Gkanogiannis <anestis@gkanogiannis.com>
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, 
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
 *
 */
package fr.cea.ig.metatarget.utils;

import java.util.concurrent.ThreadFactory;
import java.util.concurrent.atomic.AtomicLong;

public class NamedThreadFactory implements ThreadFactory {
	private static final AtomicLong THREAD_POOL_NUM = new AtomicLong(0);
	private final AtomicLong mThreadNum = new AtomicLong(0);
	private final String mPrefix;
	private final boolean mIsDaemon;
	private final long mPoolNum;

	public NamedThreadFactory(String pPrefix) {
		this(pPrefix, true);
	}

	public NamedThreadFactory(String pPrefix, boolean pIsDaemon) {
		mIsDaemon = pIsDaemon;
		mPrefix = pPrefix;
		mPoolNum = THREAD_POOL_NUM.incrementAndGet();
	}

	@Override
	public Thread newThread(Runnable r) {
		Thread t = new Thread(r, mPrefix + "-" + mPoolNum + "-Thread-" + mThreadNum.incrementAndGet());
		if (t.isDaemon() != mIsDaemon)
			t.setDaemon(mIsDaemon);

		return t;
	}
}
